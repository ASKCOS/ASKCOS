# Deploy website - first copy correct settings over
cp askcos_site/settings_deployment.py askcos_site/settings.py
if ! pidof apache2 > /dev/null
then
    python manage.py runserver
else
    sudo service apache2 reload
fi