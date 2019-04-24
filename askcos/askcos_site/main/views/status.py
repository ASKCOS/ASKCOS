from django.shortcuts import render

def status(request):
    context = {}
    return render(request, 'status.html', context)