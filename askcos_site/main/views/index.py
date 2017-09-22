from django.shortcuts import render, HttpResponse, redirect
from django.contrib.auth.decorators import login_required

@login_required
def index(request, err=None):
    '''
    Homepage
    '''
    print('{} loaded the index page!'.format(request.user))
    return render(request, 'index.html', {'err': err})