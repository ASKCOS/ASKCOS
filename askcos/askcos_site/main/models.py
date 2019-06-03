from django.db import models
from django.contrib.auth.models import User

class SavedResults(models.Model):
    user = models.ForeignKey(User)
    description = models.CharField(max_length=1000)
    created = models.DateTimeField('created on')
    dt = models.CharField(max_length=200)
    fpath = models.CharField(max_length=500)

class BlacklistedReactions(models.Model):
    user = models.ForeignKey(User)
    description = models.CharField(max_length=1000)
    created = models.DateTimeField('created on')
    dt = models.CharField(max_length=200)
    smiles = models.CharField(max_length=5000)
    active = models.BooleanField()

class BlacklistedChemicals(models.Model):
    user = models.ForeignKey(User)
    description = models.CharField(max_length=1000)
    created = models.DateTimeField('created on')
    dt = models.CharField(max_length=200)
    smiles = models.CharField(max_length=5000)
    active = models.BooleanField()
