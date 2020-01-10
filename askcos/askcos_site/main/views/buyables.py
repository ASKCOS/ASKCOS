from django.shortcuts import render
from askcos_site.main.views.users import can_modify_buyables

def buyables(request):
    return render(request, 'buyables.html', {'can_modify_buyables': can_modify_buyables(request)})
