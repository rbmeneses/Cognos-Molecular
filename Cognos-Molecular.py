import streamlit as st
import cirpy
import rdkit.Chem as Chem
from rdkit.Chem import Draw, AllChem
import io
import google.generativeai as genai
import re
import sqlite3
import pandas as pd
from datetime import datetime

# --- Configuração da Página ---
st.set_page_config(page_title="Cognos Molecular AI - V2", page_icon="🧬", layout="wide")

# --- Gerenciador de Banco de Dados (SQLite) ---
class DatabaseManager:
    def __init__(self, db_name="moleculas.db"):
        self.conn = sqlite3.connect(db_name, check_same_thread=False)
        self.create_table()

    def create_table(self):
        cursor = self.conn.cursor()
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS library (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT,
                smiles TEXT,
                formula TEXT,
                created_at TEXT
            )
        ''')
        self.conn.commit()

    def save_molecule(self, name, smiles, formula):
        try:
            cursor = self.conn.cursor()
            date_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            cursor.execute('''
                INSERT INTO library (name, smiles, formula, created_at)
                VALUES (?, ?, ?, ?)
            ''', (name, smiles, formula, date_now))
            self.conn.commit()
            return True
        except Exception as e:
            st.error(f"Erro ao salvar no BD: {e}")
            return False

    def get_all_molecules(self):
        return pd.read_sql_query("SELECT * FROM library ORDER BY id DESC", self.conn)

# --- Inteligência Artificial (Gemini) ---
class GeminiService:
    def __init__(self, api_key):
        self.api_key = api_key
        self.enabled = bool(api_key)
        if self.enabled:
            genai.configure(api_key=api_key)
            self.model = genai.GenerativeModel('gemini-2.5-flash')

    def get_smiles_from_text(self, text):
        if not self.enabled: return None
        try:
            prompt = (
                f"Tarefa: Converter texto em SMILES. Entrada: '{text}'. "
                f"Retorne APENAS a string SMILES. Se não souber, retorne INVALID."
            )
            response = self.model.generate_content(prompt)
            return response.text.strip().replace("`", "")
        except: return None

    def generate_name_for_structure(self, smiles):
        """
        Analisa a estrutura. Se existir, dá o nome real.
        Se for inventada, cria um nome IUPAC ou descritivo plausível.
        """
        if not self.enabled: return "Estrutura Desconhecida (Sem IA)"
        try:
            prompt = (
                f"Atue como um químico sênior. Analise este código SMILES: '{smiles}'. "
                f"1. Se for uma molécula real conhecida, qual o nome comum (em Português)? "
                f"2. Se for uma estrutura teórica ou com elementos curinga (*, M, R), invente um nome descritivo científico plausível baseado na estrutura. "
                f"Responda APENAS com o nome. Nada mais."
            )
            response = self.model.generate_content(prompt)
            return response.text.strip().replace("**", "")
        except Exception as e:
            return "Erro na Nomenclatura IA"

# --- Fábrica de Moléculas (Modo Nuclear/Criativo) ---
class MoleculeFactory:
    @staticmethod
    def create_mol(input_str, is_smiles=False):
        mol = None
        if not is_smiles: return None

        strategies = [
            lambda s: Chem.MolFromSmiles(s),
            lambda s: Chem.MolFromSmiles(s, sanitize=False),
            lambda s: MoleculeFactory._clean_parse(s, fix_elements=True),
            lambda s: MoleculeFactory._clean_parse(s, fix_elements=True, fix_aromaticity=True),
            lambda s: MoleculeFactory._clean_parse(s, fix_elements=True, fix_aromaticity=True, remove_rings=True)
        ]

        for strategy in strategies:
            try:
                mol = strategy(input_str)
                if mol:
                    try: mol.UpdatePropertyCache(strict=False)
                    except: pass
                    break
            except: continue
        return mol

    @staticmethod
    def _clean_parse(smiles, fix_elements=False, fix_aromaticity=False, remove_rings=False):
        s = smiles
        if remove_rings: s = re.sub(r'\d+', '', s)
        if fix_aromaticity:
            s = s.upper()
            s = s.replace("CL", "Cl").replace("BR", "Br").replace("NA", "Na")
        if fix_elements:
            for bad in ['M', 'R', 'X', 'Z', 'J', 'Q', 'A', 'D', 'E', 'G', 'L', 'T']:
                s = re.sub(fr'\b{bad}\b', '*', s)
                s = s.replace(f"={bad}", "=*").replace(f"#{bad}", "#*").replace(f"({bad})", "(*)")
        return Chem.MolFromSmiles(s, sanitize=False)

# --- Resolver Principal ---
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
        elif input_type == 'CAS Number':
            try: smiles = cirpy.resolve(value, 'smiles')
            except: pass

        if not smiles or "INVALID" in smiles: return None, None, None
        
        mol = MoleculeFactory.create_mol(smiles, is_smiles=True)
        
        # Gera nome sugerido
        suggested_name = "Desconhecido"
        if mol:
            # Tenta Cirpy reverso primeiro (Internet Database)
            try: 
                suggested_name = cirpy.resolve(smiles, 'names')[0]
            except:
                # Se falhar, usa a IA para inventar/buscar
                suggested_name = self.gemini.generate_name_for_structure(smiles)

        return mol, smiles, suggested_name

# --- Interface Gráfica ---
def main():
    st.title("🧪 Cognos Molecular AI - Database Edition")
    
    # Inicializa BD
    db = DatabaseManager()

    with st.sidebar:
        st.header("⚙️ Configurações")
        api_key = st.text_input("Gemini API Key", type="password")
        
        st.divider()
        st.header("📚 Biblioteca Salva")
        df = db.get_all_molecules()
        if not df.empty:
            st.dataframe(df[['name', 'formula', 'created_at']], hide_index=True)
            
            # Botão para baixar CSV
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Baixar Planilha Completa", csv, "biblioteca_molecular.csv", "text/csv")
        else:
            st.info("Nenhuma molécula salva ainda.")

    # Área Principal
    col1, col2 = st.columns([1, 4])
    with col1:
        input_type = st.selectbox("Tipo de Entrada", ["Nome / Descrição", "SMILES (Raw)", "CAS Number"])
    with col2:
        ph = "Ex: Sem api usar somente nome oficial em ingles BETA-NICOTINAMIDE ADENINE DINUCLEOTIDE com api beta-nicotinamida adenina dinucleotídeo ou gluconolactona  ou paracetamol" if input_type == "Nome / Descrição" else "Ex: CC(=M)Nc7ccc(cc5)O"
        user_input = st.text_input("Entrada", placeholder=ph)

    if st.button("Analisar & Gerar 🧬", type="primary", use_container_width=True):
        if not user_input: return

        gemini_svc = GeminiService(api_key)
        resolver = ChemicalResolver(gemini_svc)

        with st.spinner("Processando estrutura e buscando identidade..."):
            mol, final_smiles, suggested_name = resolver.resolve(input_type, user_input)

            if mol:
                st.session_state['current_mol'] = mol
                st.session_state['current_smiles'] = final_smiles
                st.session_state['suggested_name'] = suggested_name
            else:
                st.error("Não foi possível processar a estrutura.")

    # Exibição e Salvamento
    if 'current_mol' in st.session_state:
        st.divider()
        mol = st.session_state['current_mol']
        final_smiles = st.session_state['current_smiles']
        
        # Cálculo da Fórmula
        try: formula_str = AllChem.CalcMolFormula(mol)
        except: formula_str = "Indefinida (Estrutura Aberta)"

        c1, c2 = st.columns([1, 1])
        
        with c1:
            img = Draw.MolToImage(mol, size=(500, 500))
            st.image(img, caption="Renderização 2D", use_container_width=True)

        with c2:
            st.subheader("📝 Dados & Salvamento")
            
            # Campo de nome editável (preenchido pela IA)
            name_val = st.session_state.get('suggested_name', 'Nova Molécula')
            final_name = st.text_input("Nome da Molécula (Sugerido pela IA):", value=name_val)
            
            st.text_area("SMILES", final_smiles, height=70)
            st.info(f"Fórmula Química: **{formula_str}**")

            col_save, col_dl = st.columns(2)
            with col_save:
                if st.button("💾 Salvar no Banco de Dados", use_container_width=True):
                    if db.save_molecule(final_name, final_smiles, formula_str):
                        st.toast(f"Salvo: {final_name}", icon="✅")
                        st.rerun() # Atualiza a sidebar

            with col_dl:
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                st.download_button("📥 Download PNG", buf.getvalue(), f"{final_name}.png", "image/png", use_container_width=True)

if __name__ == "__main__":

    main()
