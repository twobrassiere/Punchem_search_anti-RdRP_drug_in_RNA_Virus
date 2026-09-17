import os
import time
import logging
import requests
import pandas as pd
from datetime import datetime
from typing import List, Dict, Set

import json
try:
    import gspread
    from google.oauth2.service_account import Credentials
    GSPREAD_AVAILABLE = True
except ImportError:
    GSPREAD_AVAILABLE = False

# 設定 Logging 紀錄，方便除錯與追蹤每日任務
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)

# 參考 AntiviralDB 收錄之完整 RNA 病毒種類清單 (包含 +ssRNA, -ssRNA, dsRNA 病毒)
ANTIVIRALDB_RNA_VIRUS_LIST = [
    "SARS-CoV-2",
    "SARS-CoV",
    "MERS-CoV",
    "Human coronavirus 229E",
    "Human coronavirus NL63",
    "Hepatitis C virus",
    "Dengue virus",
    "Zika virus",
    "West Nile virus",
    "Yellow fever virus",
    "Japanese encephalitis virus",
    "Enterovirus 71",
    "Poliovirus",
    "Rhinovirus",
    "Hepatitis A virus",
    "Coxsackievirus",
    "Norovirus",
    "Chikungunya virus",
    "Rubella virus",
    "Influenza A virus",
    "Influenza B virus",
    "Respiratory syncytial virus",
    "Nipah virus",
    "Measles virus",
    "Mumps virus",
    "Rabies virus",
    "Ebola virus",
    "Marburg virus",
    "Lassa virus",
    "Crimean-Congo hemorrhagic fever virus",
    "Rotavirus"
]

PUBCHEM_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
NCBI_ESEARCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def search_cids_via_esearch(query: str, retmax: int = 200) -> List[int]:
    """
    使用 NCBI Entrez ESearch API 進行靈活的文字關鍵字搜尋
    """
    params = {
        "db": "pccompound",
        "term": query,
        "retmode": "json",
        "retmax": retmax
    }
    
    try:
        response = requests.get(NCBI_ESEARCH_BASE, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            return [int(cid) for cid in id_list]
    except Exception as e:
        logging.error(f"ESearch 失敗 [{query}]: {e}")
        
    return []

def search_cids_via_pug_name(term: str) -> List[int]:
    """
    使用 PubChem PUG REST 按名稱搜尋 CID
    """
    url = f"{PUBCHEM_PUG_BASE}/compound/name/{requests.utils.quote(term)}/cids/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return data.get("IdentifierList", {}).get("CID", [])
    except Exception as e:
        logging.debug(f"PUG Name 搜尋無結果或失敗 [{term}]: {e}")
        
    return []


def fetch_compound_properties(cids: List[int]) -> List[Dict]:
    """
    分批從 PubChem 抓取化合物的詳細化學性質 (SMILES, MW, Formula, Title 等)
    遵守 PubChem Rate Limit Rules (每秒不超過 5 次請求)
    """
    if not cids:
        return []
    
    properties_to_fetch = "Title,IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChIKey"
    chunk_size = 100  # 避免 URL 長度過長
    all_properties = []
    
    for i in range(0, len(cids), chunk_size):
        sub_cids = ",".join(map(str, cids[i:i + chunk_size]))
        url = f"{PUBCHEM_PUG_BASE}/compound/cid/{sub_cids}/property/{properties_to_fetch}/JSON"
        
        try:
            res = requests.get(url, timeout=15)
            if res.status_code == 200:
                props = res.json().get("PropertyTable", {}).get("Properties", [])
                all_properties.extend(props)
            else:
                logging.warning(f"無法取得 CID {sub_cids[:20]}... 屬性，狀態碼: {res.status_code}")
        except Exception as e:
            logging.error(f"抓取屬性時發生錯誤: {e}")
            
        # 遵循 Rate limit: 請求之間適當延遲
        time.sleep(0.25)
        
    return all_properties


def update_google_sheet_by_virus(virus_df_dict: Dict[str, pd.DataFrame], sheet_name: str = "RNA_Virus_RdRP_Inhibitors") -> None:
    """
    將各病毒的搜尋結果分別寫入至 Google Sheet 中，並以各病毒名稱做為工作表 (Worksheet Tab) 命名。
    """
    if not GSPREAD_AVAILABLE:
        logging.warning("未安裝 gspread / google-auth 套件。請先執行 'pip install gspread google-auth'")
        return

    creds_env = os.environ.get("GOOGLE_SERVICE_ACCOUNT_JSON")
    creds_filename = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS", "credentials.json")

    scopes = [
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive"
    ]

    try:
        # 判斷是讀取 CI/CD 環境變數還是本機 JSON 憑證檔
        if creds_env:
            creds_dict = json.loads(creds_env)
            creds = Credentials.from_service_account_info(creds_dict, scopes=scopes)
        elif os.path.exists(creds_filename):
            creds = Credentials.from_service_account_file(creds_filename, scopes=scopes)
        else:
            logging.warning("找不到 Google Service Account 金鑰，跳過 Google Sheet 上傳 (請設定 credentials.json 或 GOOGLE_SERVICE_ACCOUNT_JSON)。")
            return

        client = gspread.authorize(creds)

        # 開啟試算表 (若不存在則自動新建)
        try:
            spreadsheet = client.open(sheet_name)
        except gspread.SpreadsheetNotFound:
            spreadsheet = client.create(sheet_name)
            logging.info(f"已在 Google Drive 建立新的 Google Sheet: {sheet_name}")

        existing_worksheets = {ws.title: ws for ws in spreadsheet.worksheets()}

        # 逐一針對每一種病毒更新對應名稱的工作表 (Worksheet)
        for virus_name, df in virus_df_dict.items():
            # Google Sheet 工作表名稱限制最多 100 個字元
            ws_title = virus_name[:100]

            # 轉換 DataFrame 格式寫入試算表
            if not df.empty:
                data_to_write = [df.columns.tolist()] + df.fillna("").astype(str).values.tolist()
            else:
                data_to_write = [[
                    "PubChem CID", "TargetViruses", "Organism", "Target Name", "UniProt ID",
                    "Compound Pref Name", "Standard Type", "Standard Relation", "Standard Value",
                    "Standard Units", "Canonical SMILES", "InChIKey", "Molecular Formula", "Molecular Weight"
                ]]

            if ws_title in existing_worksheets:
                worksheet = existing_worksheets[ws_title]
                try:
                    existing_data = worksheet.get_all_values()
                except Exception as e:
                    logging.warning(f"讀取工作表 [{ws_title}] 內容失敗: {e}")
                    existing_data = []

                # 比對既有內容與最新資料是否相同
                if existing_data == data_to_write:
                    logging.info(f"工作表 [{ws_title}] 內容無變動，跳過更新")
                    continue
                else:
                    worksheet.clear()
                    worksheet.update(values=data_to_write, range_name="A1")
                    logging.info(f"工作表 [{ws_title}] 檢測到數據變更，已更新為最新版本（共 {len(df)} 筆資料）")
            else:
                rows_count = max(100, len(df) + 10)
                worksheet = spreadsheet.add_worksheet(title=ws_title, rows=str(rows_count), cols="20")
                existing_worksheets[ws_title] = worksheet
                worksheet.update(values=data_to_write, range_name="A1")
                logging.info(f"建立新工作表 [{ws_title}] 並寫入最新資料（共 {len(df)} 筆資料）")

        # 若自動產生的預設空白 Sheet1 存在且已有其他病毒工作表，則進行清理
        if "Sheet1" in existing_worksheets and len(spreadsheet.worksheets()) > 1:
            try:
                spreadsheet.del_worksheet(existing_worksheets["Sheet1"])
            except Exception:
                pass

        logging.info(f"成功將所有病毒分頁資料同步至 Google Sheet: '{sheet_name}'")
        logging.info(f"Sheet 連結: {spreadsheet.url}")

    except Exception as e:
        logging.error(f"寫入 Google Sheet 時發生錯誤: {e}")


def collect_antiviraldb_rdrp_inhibitors() -> tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    依據 AntiviralDB RNA 病毒清單，逐一對每一種病毒搜尋 PubChem RdRP 抑制劑資料並依病毒分類
    """
    virus_to_cids: Dict[str, Set[int]] = {}
    all_unique_cids: Set[int] = set()

    logging.info(f"=== 開始依據 AntiviralDB 每日搜尋 RNA 病毒 RdRP 抑制劑（共 {len(ANTIVIRALDB_RNA_VIRUS_LIST)} 種病毒）===")
    
    # 逐一對 AntiviralDB 清單中的每一種 RNA 病毒進行檢索
    for idx, virus in enumerate(ANTIVIRALDB_RNA_VIRUS_LIST, start=1):
        query_terms = [
            f"{virus} RdRP inhibitor",
            f"{virus} RNA-dependent RNA polymerase inhibitor"
        ]
        
        logging.info(f"[{idx}/{len(ANTIVIRALDB_RNA_VIRUS_LIST)}] 搜尋病毒: '{virus}'")
        
        found_cids_for_virus: Set[int] = set()
        for query in query_terms:
            cids = search_cids_via_esearch(query)
            if not cids:
                cids = search_cids_via_pug_name(query)
            
            found_cids_for_virus.update(cids)
            time.sleep(0.25)  # 遵守 PubChem API 限速規範
            
        virus_to_cids[virus] = found_cids_for_virus
        all_unique_cids.update(found_cids_for_virus)

        logging.info(f" -> 病毒 '{virus}' 共找到 {len(found_cids_for_virus)} 個相關 PubChem CID")

    logging.info(f"=== 全部病毒檢索完畢！累計 {len(all_unique_cids)} 個不重複化合物 CID ===")

    column_order = [
        "PubChem CID",
        "TargetViruses",
        "Organism",
        "Target Name",
        "UniProt ID",
        "Compound Pref Name",
        "Standard Type",
        "Standard Relation",
        "Standard Value",
        "Standard Units",
        "Canonical SMILES",
        "InChIKey",
        "Molecular Formula",
        "Molecular Weight"
    ]

    if not all_unique_cids:
        empty_df = pd.DataFrame(columns=column_order)
        return empty_df, {virus: empty_df for virus in ANTIVIRALDB_RNA_VIRUS_LIST}

    # 下載 PubChem 化合物詳細化學結構與資訊
    logging.info("開始批量抓取 PubChem 化合物詳細結構資訊...")
    compounds_data = fetch_compound_properties(list(all_unique_cids))
    pubchem_dict = {comp.get("CID"): comp for comp in compounds_data if comp.get("CID")}

    virus_df_dict: Dict[str, pd.DataFrame] = {}
    all_rows = []

    # 針對每一種病毒獨立產生其專屬的 DataFrame
    for virus in ANTIVIRALDB_RNA_VIRUS_LIST:
        cids_for_this_virus = virus_to_cids.get(virus, set())
        rows = []
        for cid in cids_for_this_virus:
            comp = pubchem_dict.get(cid, {})
            comp_info = {
                "PubChem CID": cid,
                "TargetViruses": virus,
                "Organism": virus,
                "Target Name": f"RdRP ({virus})",
                "UniProt ID": "N/A (PubChem Search)",
                "Compound Pref Name": comp.get("Title", ""),
                "Standard Type": "Keyword Match",
                "Standard Relation": "=",
                "Standard Value": "N/A",
                "Standard Units": "nM",
                "Canonical SMILES": comp.get("CanonicalSMILES", ""),
                "InChIKey": comp.get("InChIKey", ""),
                "Molecular Formula": comp.get("MolecularFormula", ""),
                "Molecular Weight": comp.get("MolecularWeight", "")
            }
            rows.append(comp_info)
            all_rows.append(comp_info)

        df_virus = pd.DataFrame(rows)
        if not df_virus.empty:
            df_virus = df_virus.reindex(columns=column_order)
        else:
            df_virus = pd.DataFrame(columns=column_order)
            
        virus_df_dict[virus] = df_virus

    df_combined = pd.DataFrame(all_rows)
    if not df_combined.empty:
        df_combined = df_combined.reindex(columns=column_order)
    else:
        df_combined = pd.DataFrame(columns=column_order)

    return df_combined, virus_df_dict


def main():
    # 執行依 AntiviralDB 清單之每日自動化搜尋
    df_combined, virus_df_dict = collect_antiviraldb_rdrp_inhibitors()
    
    if df_combined.empty:
        logging.warning("未抓取到任何資料。")
        return

    # 建立輸出資料夾與 CSV 檔
    output_dir = "data"
    os.makedirs(output_dir, exist_ok=True)
    
    today = datetime.now().strftime("%Y-%m-%d")
    csv_filename = os.path.join(output_dir, f"rdrp_inhibitors_antiviralDB_{today}.csv")
    
    # 寫入 CSV 檔案（總合資料）
    df_combined.to_csv(csv_filename, index=False, encoding="utf-8-sig")
    logging.info(f"成功將 {len(df_combined)} 筆 RdRP 抑制劑總合資料寫入至: {csv_filename}")

    # 自動覆寫同步至 Google Sheets，每個病毒獨立一個工作表 (Worksheet Tab)
    update_google_sheet_by_virus(virus_df_dict, sheet_name="RNA_Virus_RdRP_Inhibitors")

if __name__ == "__main__":
    main()