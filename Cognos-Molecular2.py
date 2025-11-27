import streamlit as st
import cirpy
import rdkit.Chem as Chem
from rdkit.Chem import Draw, AllChem
import io
import google.generativeai as genai
import re

# --- Configuração da Página ---
st.set_page_config(page_title="Cognos Molecular", page_icon="🧪", layout="wide")

# --- Classes de Serviço ---

class GeminiService:
    def __init__(self, api_key):
        self.api_key = api_key
        self.enabled = bool(api_key)
        if self.enabled:
            genai.configure(api_key=api_key)

    def get_smiles_from_text(self, text):
        if not self.enabled: 
            return None
        try:
            model = genai.GenerativeModel('gemini-2.5-flash')
            prompt = (
                f"Atue como um químico especialista. Entrada: '{text}'. "
                f"Tarefa: Identifique a molécula (aceite nomes em Português-BR, Inglês, IUPAC, marcas comerciais). "
                f"Saída: Retorne APENAS a string SMILES oficial. Sem explicações, sem markdown, sem aspas. "
                f"Se a entrada for um SMILES inválido quimicamente mas sintaticamente correto, retorne ele mesmo. "
                f"Se não conseguir identificar, retorne INVALID."
            )
            response = model.generate_content(prompt)
            result = response.text.strip().replace("`", "").replace("smiles", "").strip()
            return result if "INVALID" not in result else None
        except Exception as e:
            st.error(f"Erro na comunicação com a IA: {e}")
            return None

class MoleculeFactory:
    """
    Fábrica de Moléculas com tratamento agressivo de erros.
    """
    @staticmethod
    def create_mol(input_str, is_smiles=False):
        mol = None
        if not is_smiles: return None

        # Lista de Estratégias (da mais correta para a mais "desesperada")
        strategies = [
            # 1. Padrão: Tenta ler SMILES correto
            lambda s: Chem.MolFromSmiles(s),
            
            # 2. Leniente: Ignora erros de valência (ex: Carbono com 5 pernas)
            lambda s: Chem.MolFromSmiles(s, sanitize=False),
            
            # 3. Criativo: Substitui letras inventadas (M, R, Z) por * (Curinga)
            lambda s: MoleculeFactory._clean_parse(s, fix_elements=True),
            
            # 4. Corretivo: Fixa elementos E remove aromaticidade (c->C) para evitar anéis quebrados
            lambda s: MoleculeFactory._clean_parse(s, fix_elements=True, fix_aromaticity=True),

            # 5. NUCLEAR: Fixa tudo E remove números (Ignora anéis que não fecham, ex: "c7" sem par)
            lambda s: MoleculeFactory._clean_parse(s, fix_elements=True, fix_aromaticity=True, remove_rings=True)
        ]

        for strategy in strategies:
            try:
                mol = strategy(input_str)
                if mol:
                    try: mol.UpdatePropertyCache(strict=False)
                    except: pass
                    break
            except:
                continue
        
        return mol

    @staticmethod
    def _clean_parse(smiles, fix_elements=False, fix_aromaticity=False, remove_rings=False):
        s = smiles
        
        # Passo 1: Remove números se solicitado (Resolve erro de anel aberto "7" ou "5")
        if remove_rings:
            s = re.sub(r'\d+', '', s)

        # Passo 2: Transforma letras minúsculas em maiúsculas (Resolve erro de aromaticidade)
        if fix_aromaticity:
            # Transforma c, n, o, s, p em maiúsculas
            # Nota: Isso quebra Cl, Br, mas o RDKit costuma lidar bem se o resto estiver ok.
            # Uma abordagem mais segura é upper() em tudo se for o modo nuclear
            s = s.upper() 
            # Correção pós-upper para Cloro (CL -> Cl) e Bromo (BR -> Br) para não virarem Carbono + Leucina fictícia
            s = s.replace("CL", "Cl").replace("BR", "Br").replace("NA", "Na")

        # Passo 3: Substitui elementos fictícios (M, J, X, Z) por *
        if fix_elements:
            # Substitui letras isoladas que não são elementos comuns orgânicos
            # Regex procura letras que NÃO sejam B, C, N, O, P, S, F, I, H...
            # Simplificação: Troca M, R, X, Z, J, Q diretos
            for bad in ['M', 'R', 'X', 'Z', 'J', 'Q', 'A', 'D', 'E', 'G', 'L', 'T']:
                # Evita substituir dentro de "Cl" ou "Br" se já tratamos
                s = re.sub(fr'\b{bad}\b', '*', s) # Letra isolada
                s = s.replace(f"={bad}", "=*").replace(f"#{bad}", "#*").replace(f"({bad})", "(*)")
        
        return Chem.MolFromSmiles(s, sanitize=False)

class ChemicalResolver:
    def __init__(self, gemini_service):
        self.gemini = gemini_service

    def resolve(self, input_type, value):
        smiles = None
        if input_type == 'SMILES (Raw)':
            smiles = value.strip()
        elif input_type == 'Nome / Descrição':
            try: smiles = cirpy.resolve(value, 'smiles')
            except: pass
            if not smiles:
                smiles = self.gemini.get_smiles_from_text(value)
                if smiles: st.toast(f"Gemini identificou: {smiles}", icon="🤖")
        elif input_type == 'CAS Number':
            try: smiles = cirpy.resolve(value, 'smiles')
            except: pass

        if not smiles: return None, None
        mol = MoleculeFactory.create_mol(smiles, is_smiles=True)
        return mol, smiles

# --- Interface Gráfica ---

def main():
    st.title("🧪 Cognos Molecular AI")
    st.markdown("""
    **Modo Universal:** Desenha estruturas reais ou hipotéticas.
    Se o código contiver erros (anéis abertos, elementos inventados), o sistema tentará corrigir automaticamente.
    """)

    with st.sidebar:
        st.header("Configurações")
        api_key = st.text_input("Gemini API Key", type="password")
        st.divider()

    col1, col2 = st.columns([1, 4])
    with col1:
        input_type = st.selectbox("Tipo de Entrada", ["Nome / Descrição", "SMILES (Raw)", "CAS Number"])
    with col2:
        ph = "Ex: Paracetamol" if input_type == "Nome / Descrição" else "Ex: CC(=M)Nc7ccc(cc5)O"
        user_input = st.text_input("Entrada", placeholder=ph)

    if st.button("Gerar Estrutura 🧬", type="primary", use_container_width=True):
        if not user_input:
            st.warning("Digite algo para começar.")
            return

        gemini_svc = GeminiService(api_key)
        resolver = ChemicalResolver(gemini_svc)

        with st.spinner("Processando..."):
            mol, final_smiles = resolver.resolve(input_type, user_input)

            if mol:
                st.session_state['mol'] = mol
                st.session_state['smiles'] = final_smiles
                st.success(f"Sucesso!")
            else:
                st.error("Erro fatal: A estrutura é impossível de desenhar mesmo com correções.")

    if 'mol' in st.session_state and st.session_state['mol']:
        st.divider()
        mol = st.session_state['mol']
        
        c1, c2 = st.columns([1, 1])
        with c1:
            try:
                img = Draw.MolToImage(mol, size=(600, 600))
                st.image(img, caption="Estrutura Visualizada", use_container_width=True)
                
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                st.download_button("📥 Baixar PNG", buf.getvalue(), "estrutura.png", "image/png", use_container_width=True)
            except:
                st.error("Erro na renderização da imagem.")

        with c2:
            st.subheader("SMILES Processado")
            st.code(st.session_state['smiles'])
            try:
                wt = AllChem.CalcExactMolWt(mol)
                st.info(f"Peso Molecular Estimado: {wt:.2f} g/mol")
            except:
                st.warning("Propriedades físicas indisponíveis (Estrutura hipotética).")

if __name__ == "__main__":
    main()