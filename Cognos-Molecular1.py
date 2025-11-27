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
    """
    Gerencia a comunicação com a IA para resolver nomes em PT-BR ou descrições complexas.
    """
    def __init__(self, api_key):
        self.api_key = api_key
        self.enabled = bool(api_key)
        if self.enabled:
            genai.configure(api_key=api_key)

    def get_smiles_from_text(self, text):
        if not self.enabled: 
            return None
        try:
            model = genai.GenerativeModel('gemini-2.5-flash') # Ou gemini-pro
            # Prompt otimizado para extrair apenas o SMILES de nomes em PT-BR
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
    Cria objetos RDKit. 
    Agora com 'Modo Nuclear' que corrige aromaticidade quebrada (c -> C) e elementos fictícios.
    """
    @staticmethod
    def create_mol(input_str, is_smiles=False):
        mol = None
        if not is_smiles: return None

        # Sequência de tentativas: da mais correta para a mais "forçada"
        attempts = [
            # 1. Padrão (Química Perfeita)
            lambda s: Chem.MolFromSmiles(s),
            
            # 2. Leniente (Ignora valência impossível)
            lambda s: Chem.MolFromSmiles(s, sanitize=False),
            
            # 3. Criativo (Troca M/R/Z por *, mas mantém aromaticidade)
            lambda s: MoleculeFactory._try_creative_parsing(s, force_dearomatization=False),
            
            # 4. NUCLEAR (Troca M/R/Z por * E converte c->C, n->N para ignorar anéis abertos)
            lambda s: MoleculeFactory._try_creative_parsing(s, force_dearomatization=True)
        ]

        for attempt in attempts:
            try:
                mol = attempt(input_str)
                if mol:
                    try: mol.UpdatePropertyCache(strict=False)
                    except: pass
                    break
            except:
                continue

        return mol

    @staticmethod
    def _try_creative_parsing(smiles, force_dearomatization=False):
        # 1. Substitui elementos fictícios (M, R, X, Z, J) por '*' (Dummy Atom)
        s = smiles
        # Substituição via Regex para letras isoladas
        s = re.sub(r'\b[MRXZJ]\b', '*', s) 
        # Substituição direta para casos "grudados" como (=M)
        for bad in ['M', 'R', 'X', 'Z', 'J']:
            s = s.replace(bad, '*')

        # 2. Se ativado, remove a aromaticidade (c -> C)
        # Isso permite desenhar cadeias que o usuário digitou como 'c' minúsculo
        # mas esqueceu de fechar o anel.
        if force_dearomatization:
            # Substitui letras minúsculas comuns por maiúsculas
            # Cuidado para não quebrar Cl, Br, etc (que usam segunda letra minúscula)
            # Substituímos apenas os átomos orgânicos comuns de anéis
            for lower, upper in [('c', 'C'), ('n', 'N'), ('o', 'O'), ('s', 'S'), ('p', 'P')]:
                s = s.replace(lower, upper)

        return Chem.MolFromSmiles(s, sanitize=False)

class ChemicalResolver:
    """Orquestra a decisão entre buscar no Cirpy (rápido) ou Gemini (inteligente)."""
    def __init__(self, gemini_service):
        self.gemini = gemini_service

    def resolve(self, input_type, value):
        smiles = None
        
        # Se o usuário escolheu digitar o SMILES diretamente
        if input_type == 'SMILES (Raw)':
            smiles = value
            # Remove espaços acidentais
            smiles = smiles.strip()

        # Se o usuário digitou um nome (ex: Metanol, Paracetamol)
        elif input_type == 'Nome / Descrição':
            # 1. Tenta Cirpy primeiro (Bases oficiais)
            try:
                smiles = cirpy.resolve(value, 'smiles')
            except:
                pass

            # 2. Se falhar, usa o Gemini (Entende Português e erros de digitação)
            if not smiles:
                smiles = self.gemini.get_smiles_from_text(value)
                if smiles: 
                    st.toast(f"Gemini identificou a estrutura!", icon="🤖")

        # Se o usuário digitou CAS
        elif input_type == 'CAS Number':
            try:
                smiles = cirpy.resolve(value, 'smiles')
            except:
                pass

        if not smiles:
            return None, None

        # Gera o objeto visual (Mol) permitindo erros
        mol = MoleculeFactory.create_mol(smiles, is_smiles=True)
        return mol, smiles

# --- Interface Gráfica (Streamlit) ---

def main():
    st.title("🧪 Cognos Molecular AI")
    st.markdown("""
    Ferramenta de visualização molecular. 
    - **Nomes em Português:** Reconhecidos via Google Gemini.
    - **SMILES Livre:** Desenha estruturas mesmo com erros de valência (ex: Carbono com 5 ligações).
    """)

    # Sidebar para API Key
    with st.sidebar:
        st.header("Configurações")
        api_key = st.text_input("Gemini API Key", type="password", help="Necessária para converter nomes em PT-BR.")
        st.caption("Sem a chave, apenas buscas exatas (inglês) e SMILES diretos funcionarão.")
        st.divider()
        st.info("Desenvolvido para análise flexível de estruturas químicas.")

    # Área Principal
    col1, col2 = st.columns([1, 4])
    with col1:
        input_type = st.selectbox("Tipo de Entrada", ["Nome / Descrição", "SMILES (Raw)", "CAS Number"])
    with col2:
        placeholder_txt = "Ex: Paracetamol, Metanol" if input_type == "Nome / Descrição" else "Ex: CC(=O)Nc1ccc(cc1)O"
        user_input = st.text_input("Entrada de Dados", placeholder=placeholder_txt)

    # Botão de Ação
    if st.button("Gerar Estrutura 🧬", type="primary", use_container_width=True):
        if not user_input:
            st.warning("Por favor, insira um nome ou código SMILES.")
            return

        gemini_svc = GeminiService(api_key)
        resolver = ChemicalResolver(gemini_svc)

        with st.spinner("Sintetizando estrutura digital..."):
            mol, final_smiles = resolver.resolve(input_type, user_input)

            if mol:
                st.session_state['mol'] = mol
                st.session_state['smiles'] = final_smiles
                st.success(f"Estrutura gerada para: **{user_input}**")
            else:
                st.error(f"Não foi possível gerar a estrutura. Se for um SMILES manual, verifique se os símbolos dos elementos existem na tabela periódica.")

    # Renderização do Resultado
    if 'mol' in st.session_state and st.session_state['mol']:
        st.divider()
        
        mol = st.session_state['mol']
        smiles_code = st.session_state['smiles']

        c1, c2 = st.columns([1, 1])

        with c1:
            # Desenha a imagem
            try:
                img = Draw.MolToImage(mol, size=(600, 600))
                st.image(img, caption="Representação 2D", use_container_width=True)
                
                # Botão de Download
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                byte_im = buf.getvalue()
                
                st.download_button(
                    label="📥 Download Imagem (PNG)",
                    data=byte_im,
                    file_name="molecula_cognos.png",
                    mime="image/png",
                    use_container_width=True
                )
            except Exception as e:
                st.error(f"Erro ao renderizar imagem: {e}")

        with c2:
            st.subheader("Dados Químicos")
            st.code(smiles_code, language="text")
            
            # Tenta calcular propriedades (pode falhar se a molécula for "impossível")
            try:
                formula = AllChem.CalcMolFormula(mol)
                weight = AllChem.CalcExactMolWt(mol)
                st.markdown(f"**Fórmula:** `{formula}`")
                st.markdown(f"**Peso Molecular:** `{weight:.3f} g/mol`")
            except:
                st.warning("⚠️ **Modo Leniente Ativo:** A estrutura contém valências ou ligações quimicamente inválidas. As propriedades físicas não podem ser calculadas, mas o desenho foi gerado conforme solicitado.")

if __name__ == "__main__":
    main()