import os

def customization(request):
    customizable = {
        'organization': os.environ.get('ORGANIZATION', 'MIT'),
        'contact_email': os.environ.get('CONTACT_EMAIL', 'ccoley@mit.edu'),
        'email_enabled': bool(os.environ.get('EMAIL_ENABLED', False)),
        'version_number': os.environ.get('VERSION_NUMBER', '0.x'),
    }
    return customizable