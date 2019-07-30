from django.shortcuts import render, HttpResponse, redirect
from django.contrib.auth.decorators import login_required

def index(request, err=None):
    '''
    Homepage
    '''
    print('{} loaded the index page!'.format(request.user))
    return render(request, 'index.html', {'err': err})

#@login_required
def modules(request, err=None):
    '''list of modules
    '''

    return render(request, 'modules.html', {'err': err})

#@login_required
def faq(request, err=None):
    return render(request, 'faq.html', {'err': err})

#@login_required
def tutorial(request, err=None):
    return render(request, 'tutorial.html', {'err': err})
