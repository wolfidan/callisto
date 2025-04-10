from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload
import os
import pickle

SCOPES = ['https://www.googleapis.com/auth/drive.file']

def authenticate():
    creds = None
    # The file token.pickle stores the user's access and refresh tokens
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    
    # If credentials are invalid or don't exist, let the user log in
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'control\credentials.json', 
                SCOPES,
                redirect_uri='http://localhost:8080'
            )
            creds = flow.run_local_server(
                port=8080,
                #access_type='offline',
                prompt='consent'
            )
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)
    
    return creds

def upload_file(filename, filepath, mimetype='image/png'):
    creds = authenticate()
    service = build('drive', 'v3', credentials=creds)
    
    file_metadata = {'name': filename}
    media = MediaFileUpload(filepath, mimetype=mimetype)
    
    try:
        file = service.files().create(
            body=file_metadata,
            media_body=media,
            fields='id'
        ).execute()
        print(f"File uploaded successfully. File ID: {file.get('id')}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    upload_file('measurement.png', 'C:\\xrt\\this_image')