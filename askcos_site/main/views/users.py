from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from pymongo.message import bson
from bson.objectid import ObjectId
from datetime import datetime
import os

from ..models import SavedResults, BlacklistedReactions, BlacklistedChemicals

can_control_robot = lambda request: request.user.get_username() in ['ccoley']

def can_view_reaxys(request):
    return request.user.is_authenticated and request.user.groups.filter(name='reaxys_view').exists()

def log_this_request(method):
    def f(*args, **kwargs):
        try:
            print('User %s requested view %s with args %r and kwargs %r' % \
                args[0].user.get_username(), method.__name__, args, kwargs)
        except Exception as e:
            print(e)
        return method(*args, **kwargs)
    return f


def login(request):
    '''
    User login
    '''
    return django.contrib.auth.views.login(request, template_name='login.html')

def logout(request):
    '''
    User logout
    '''
    return django.contrib.auth.views.logout(request, template_name='logout.html')

@login_required
def user_saved_results(request, err=None):
    saved_results = SavedResults.objects.filter(user=request.user)
    return render(request, 'saved_results.html', {'saved_results':saved_results, 'err': err})

@login_required
def user_saved_results_id(request, _id=-1):
    saved_result = SavedResults.objects.filter(user=request.user, id=_id)
    if saved_result.count() == 0:
        return user_saved_results(request, err='Could not find that ID')
    with open(saved_result[0].fpath, 'r') as fid:
        html = fid.read().decode('utf8')
    return render(request, 'saved_results_id.html', 
        {'saved_result':saved_result[0], 'html':html, 'dt': saved_result[0].dt})

@login_required
def user_saved_results_del(request, _id=-1):
    obj = SavedResults.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        os.remove(obj[0].fpath)
        obj[0].delete()
    return user_saved_results(request, err=None)

@login_required
def ajax_user_save_page(request):
    html = request.POST.get('html', None)
    desc = request.POST.get('desc', 'no description')
    dt   = request.POST.get('datetime', datetime.utcnow().strftime('%B %d, %Y %H:%M:%S %p UTC'))
    if html is None:
        print('Got None html')
        data = {'err': 'Could not get HTML to save'}
        return JsonResponse(data)
    print('Got request to save a page')
    now = datetime.now()
    unique_str = '%i.txt' % hash((now, request.user))
    fpath = os.path.join(settings.LOCAL_STORAGE['user_saves'], unique_str)
    try:
        with open(fpath, 'w') as fid:
            fid.write(html.encode('utf8'))
        print('Wrote to {}'.format(fpath))
        obj = SavedResults.objects.create(user=request.user, 
            description=desc,
            dt=dt,
            created=now,
            fpath=fpath)
        print('Created saved object {}'.format(obj.id))
    except Exception as e:
        print(e)
        print('Could not write to file {}?'.format(fpath))
        return JsonResponse({'err': 'Could not write to file {}?'.format(fpath)})

    return JsonResponse({'err': False})

@login_required
def user_blacklisted_reactions(request, err=None):
    blacklisted_reactions = BlacklistedReactions.objects.filter(user=request.user)
    return render(request, 'blacklisted_reactions.html', {'blacklisted_reactions':blacklisted_reactions, 'err': err})

@login_required
def user_blacklisted_reactions_del(request, _id=-1):
    obj = BlacklistedReactions.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        obj[0].delete()
    return user_blacklisted_reactions(request, err=None)

@login_required
def ajax_user_blacklist_reaction(request):
    smiles = request.POST.get('smiles', None)
    desc = request.POST.get('desc', 'no description')
    dt   = request.POST.get('datetime', datetime.utcnow().strftime('%B %d, %Y %H:%M:%S %p UTC'))
    if smiles is None:
        print('Got None smiles')
        data = {'err': 'Could not get reaction SMILES to save'}
        return JsonResponse(data)
    print('Got request to block a reaction')
    obj = BlacklistedReactions.objects.create(user=request.user, 
        description=desc,
        dt=dt,
        created=datetime.now(),
        smiles=smiles,
        active=True)
    print('Created blacklisted reaction object {}'.format(obj.id))
    return JsonResponse({'err': False})

@login_required 
def ajax_user_activate_reaction(request):
    _id = request.GET.get('id', -1)
    obj = BlacklistedReactions.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        obj[0].active = True 
        obj[0].save()
        return JsonResponse({'err': False})
    return JsonResponse({'err': 'Could not activate?'})

@login_required 
def ajax_user_deactivate_reaction(request):
    _id = request.GET.get('id', -1)
    obj = BlacklistedReactions.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        obj[0].active = False 
        obj[0].save()
        return JsonResponse({'err': False})
    return JsonResponse({'err': 'Could not deactivate?'})




@login_required
def user_blacklisted_chemicals(request, err=None):
    blacklisted_chemicals = BlacklistedChemicals.objects.filter(user=request.user)
    return render(request, 'blacklisted_chemicals.html', {'blacklisted_chemicals':blacklisted_chemicals, 'err': err})

@login_required
def user_blacklisted_chemicals_del(request, _id=-1):
    obj = BlacklistedChemicals.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        obj[0].delete()
    return user_blacklisted_chemicals(request, err=None)

@login_required
def ajax_user_blacklist_chemical(request):
    smiles = request.POST.get('smiles', None)
    desc = request.POST.get('desc', 'no description')
    dt   = request.POST.get('datetime', datetime.utcnow().strftime('%B %d, %Y %H:%M:%S %p UTC'))
    if smiles is None:
        print('Got None smiles')
        data = {'err': 'Could not get chemical SMILES to save'}
        return JsonResponse(data)
    print('Got request to block a chemical')
    obj = BlacklistedChemicals.objects.create(user=request.user, 
        description=desc,
        dt=dt,
        created=datetime.now(),
        smiles=smiles,
        active=True)
    print('Created blacklisted chemical object {}'.format(obj.id))
    return JsonResponse({'err': False})

@login_required 
def ajax_user_activate_chemical(request):
    _id = request.GET.get('id', -1)
    obj = BlacklistedChemicals.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        obj[0].active = True 
        obj[0].save()
        return JsonResponse({'err': False})
    return JsonResponse({'err': 'Could not activate?'})

@login_required 
def ajax_user_deactivate_chemical(request):
    _id = request.GET.get('id', -1)
    obj = BlacklistedChemicals.objects.filter(user=request.user, id=_id)
    if len(obj) == 1:
        obj[0].active = False 
        obj[0].save()
        return JsonResponse({'err': False})
    return JsonResponse({'err': 'Could not deactivate?'})