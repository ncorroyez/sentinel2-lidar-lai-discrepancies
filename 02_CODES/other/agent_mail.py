import sys
# Exclude system packages folder to avoid conflicts
if "/usr/lib/python3/dist-packages" in sys.path:
    sys.path.remove("/usr/lib/python3/dist-packages")

import imaplib
import email
from email.header import decode_header
from email.utils import parsedate_to_datetime
import google.generativeai as genai
from bs4 import BeautifulSoup
import datetime
import os
import time
from pathlib import Path
from collections import defaultdict

# ==========================================
# ⚙️ CONFIGURATION
# ==========================================

# Insert your API key here
GOOGLE_API_KEY = "AIzaSyBsuhZjtJRSe2TxgejeWSVIL34tW-FuL0M"

# Your email accounts
ACCOUNTS = [
    {
        "name": "Gmail Nathan",
        "server": "imap.gmail.com",
        "user": "nathan.corroyez14@gmail.com",
        "pass": "pudg eedn nurl dbtm" 
    }
]

# Timeframe to process (in days)
DAYS_TO_CHECK = 7

# Output directory
base_dir = Path("/home/corroyez/Documents/Obsidian/LikhoVault/00_System/02_Introspection/02.2_Recaps")

# ==========================================
# 🛠️ FUNCTIONS
# ==========================================

def clean_body(html_content):
    """Cleans HTML but keeps structure for Gemini Pro."""
    try:
        soup = BeautifulSoup(html_content, "html.parser")
        # Remove styles and scripts to save tokens
        for script in soup(["script", "style"]):
            script.extract()
        text = soup.get_text(separator=' ', strip=True)
        return text[:50000]
    except:
        return ""
    
def get_imap_date_string(days_ago):
    """Calculates the 'SINCE' date in IMAP format (DD-Mon-YYYY)."""
    english_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    cutoff_date = datetime.date.today() - datetime.timedelta(days=days_ago)
    return f"{cutoff_date.day:02d}-{english_months[cutoff_date.month - 1]}-{cutoff_date.year}"

def get_emails_grouped_by_date(account):
    """Fetches unseen emails and groups them by their parsed date."""
    print(f"📡 Connecting to {account['name']}...")
    grouped_emails = defaultdict(list)
    
    try:
        mail = imaplib.IMAP4_SSL(account['server'])
        mail.login(account['user'], account['pass'])
        mail.select("inbox")
        
        # Fetch UNSEEN emails
        imap_search_date = get_imap_date_string(DAYS_TO_CHECK)
        search_query = f'(UNSEEN SINCE "{imap_search_date}")'
        status, messages = mail.search(None, search_query)
        email_ids = messages[0].split()
        
        if not email_ids:
            print(f"   -> Nothing new on {account['name']}.")
            return grouped_emails

        print(f"   -> {len(email_ids)} emails fetched.")

        for e_id in email_ids:
            res, msg_data = mail.fetch(e_id, '(RFC822)')
            for response_part in msg_data:
                if isinstance(response_part, tuple):
                    msg = email.message_from_bytes(response_part[1])
                    
                    # Subject
                    subject, encoding = decode_header(msg["Subject"])[0]
                    if isinstance(subject, bytes):
                        subject = subject.decode(encoding if encoding else "utf-8")
                    
                    # Sender
                    from_ = msg.get("From")

                    # Parse exact Date
                    date_str = msg.get("Date")
                    try:
                        parsed_date = parsedate_to_datetime(date_str).date()
                    except Exception:
                        parsed_date = datetime.date.today()

                    # Body
                    body = ""
                    if msg.is_multipart():
                        for part in msg.walk():
                            if part.get_content_type() == "text/plain":
                                body += part.get_payload(decode=True).decode(errors='ignore')
                            elif part.get_content_type() == "text/html":
                                body += clean_body(part.get_payload(decode=True).decode(errors='ignore'))
                    else:
                        body = msg.get_payload(decode=True).decode(errors='ignore')
                    
                    extracted_text = f"\n--- EMAIL ({account['name']}) ---\nDATE: {date_str}\nDE: {from_}\nSUJET: {subject}\nCONTENU: {body}\n"
                    grouped_emails[parsed_date].append(extracted_text)
        
        mail.close()
        mail.logout()
        return grouped_emails

    except Exception as e:
        print(f"❌ Error on {account['name']}: {e}")
        return grouped_emails

def process_chunk(chunk_emails, model):
    """Generates an intermediate detailed analysis."""
    prompt_chunk = f"""
    Tu es un analyste de veille stratégique.
    Voici un lot d'emails bruts.
    
    TACHE : Catégorise les mails
    1. Pour les emails **PERSONNEL, URGENT** : Résume en 1 phrase ou 2.
    2. Pour les emails **NEWSLETTERS, ARTICLES, VEILLE, FWD, WRITINGS, BLOGS** : C'est le plus important. 
       - Ne fais PAS de résumé générique.
       - EXTRAIS les titres des sections, les liens importants, et les idées clés sous forme de liste à puces.
       - Je veux que tu détailles au possible.
    3. Pour les emails **SHOPPING, PUBS, SPAM**: Résume en 1 phrase courte.
       
    ---
    EMAILS BRUTS :
    {chunk_emails}
    """
    try:
        response = model.generate_content(prompt_chunk)
        return response.text
    except Exception as e:
        print(f"❌ Chunk error: {e}")
        return ""

# ==========================================
# 🧠 MAIN
# ==========================================

def main():
    # 1. Fetch and group all emails by date
    all_emails_by_date = defaultdict(list)
    for acc in ACCOUNTS:
        acc_emails = get_emails_grouped_by_date(acc)
        for date_key, email_list in acc_emails.items():
            all_emails_by_date[date_key].extend(email_list)
    
    if not all_emails_by_date:
        print("✅ No emails to process.")
        return

    # 2. IA Configuration
    print("🤖 Initializing Gemini 2.5 Flash...")
    genai.configure(api_key=GOOGLE_API_KEY)
    model = genai.GenerativeModel('gemini-2.5-flash')
    
    # 3. Process day by day
    # Sort dates to process oldest first (or newest first)
    sorted_dates = sorted(all_emails_by_date.keys())
    
    for processing_date in sorted_dates:
        print(f"\n=============================================")
        print(f"📅 PROCESSING DATE: {processing_date}")
        print(f"=============================================")
        
        emails_list = all_emails_by_date[processing_date]
        CHUNK_SIZE = 10
        intermediate_summaries = []
        total_chunks = (len(emails_list) + CHUNK_SIZE - 1) // CHUNK_SIZE
        
        # 4. Chunk processing for the specific date
        for i in range(0, len(emails_list), CHUNK_SIZE):
            chunk = emails_list[i:i + CHUNK_SIZE]
            chunk_dump = "\n".join(chunk) 
            
            print(f"   -> Detailed analysis for chunk {i//CHUNK_SIZE + 1} / {total_chunks}...")
            
            summary = process_chunk(chunk_dump, model)
            intermediate_summaries.append(summary)
            
            # Security pause for API quotas
            time.sleep(120) 
            
        final_input = "\n\n".join(intermediate_summaries)

        # 5. Final Synthesis for the specific date
        print(f"🤖 Generating Final Report for {processing_date}...")
        
        final_prompt = f"""
        Tu es un Consultant en Veille Stratégique.
        
        Voici tes notes d'analyse intermédiaires pour la journée du {processing_date}:
        {final_input}

        Génère un rapport MARKDOWN strict.
        
        CONSIGNES DE TRI :
        - Dans chaque catégorie, **DETAILLE chaque email puis GROUPE les explications par EXPEDITEUR**.
        - Le rapport concerne uniquement les mails du {processing_date}, tu n'as pas besoin de répéter la date à chaque ligne, concentre-toi sur le contenu.

        STRUCTURE À RESPECTER :

        # 🗞️ Daily Digest & Veille du {processing_date}

        ## 🚨 Urgences & Personnel (Action Requise)
        * [ ] **[Nom Expéditeur]** : Sujet - Action à faire.

        ## 📰 VEILLE & NEWSLETTERS (DÉTAILLÉ)
        *Groupe ici par Newsletter/Expéditeur si plusieurs mails du même.*
        
        ### 👤 [Nom de l'Expéditeur / Newsletter]
        * **Sujet principal :** ...
        * **Contenu détaillé :**
            * *Article / Point 1 :* Résumé détaillé.
            * *Article / Point 2 :* Résumé détaillé.
            * *Insight :* Chiffre ou info clé.

        ## 📉 Administratif & Divers
        * **[Expéditeur]** : Résumé bref.

        ## 🗑️ Spams & Pubs ignorées
        * Liste rapide des expéditeurs.
        """
        
        try:
            response = model.generate_content(final_prompt)
            
            # Create specific file for this date
            output_file = base_dir / f"Resume_Mails_{processing_date}.md"
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(response.text)
            
            print(f"✨ Done! Detailed report available at: {output_file}")
            
        except Exception as e:
            print(f"❌ Final error for {processing_date}: {e}")

if __name__ == "__main__":
    main()