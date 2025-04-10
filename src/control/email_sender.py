import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from datetime import datetime

# Email configuration
sender = "callisto.meteoswiss@outlook.com"
recipients = ["andrea.francesco.battaglia@irsol.usi.ch", "afbattaglia@proton.me"]
smtp_server = "smtp-mail.outlook.com"
smtp_port = 587
password = "t34min$I11"  # Store this securely!

# Create the email
def send_email():
    # Create message container
    msg = MIMEMultipart()
    msg['From'] = sender
    msg['To'] = ", ".join(recipients)
    msg['Subject'] = f"[CALLISTO@MeteoSwiss] Daily solar flux measurements: {datetime.now().strftime('%Y-%m-%d')}"

    # Add body
    body = "Please find attached today's updated solar flux summary plot, which includes the latest measurements. \n\nBest,\nCALLISTO"
    msg.attach(MIMEText(body, 'plain'))

    # Add attachment
    #with open(r'c:\xrt\output\data\proc_daily_flux', 'rb') as f:
    #    img = MIMEImage(f.read())
    #    img.add_header('Content-Disposition', 'attachment', filename='daily_fluxes_month.png')
    #    msg.attach(img)

    # Send the email
    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.starttls()
        server.login(sender, password)
        server.send_message(msg)
        print("Email sent successfully!")
    except Exception as e:
        print(f"Error sending email: {str(e)}")
    finally:
        server.quit()

if __name__ == "__main__":
    send_email()