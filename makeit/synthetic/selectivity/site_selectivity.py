from makeit.synthetic.selectivity.multitask_model import tf_predictor
import rdkit.Chem as Chem 
import sys
import os 


class Site_Predictor():
	def __init__(self):
		self.site_predictor = tf_predictor()
		self.site_predictor.build()
		print('Loaded recommendation model')
		print('### RECOMMENDER STARTED UP ###')

	def predict(self, smi):
		res = self.site_predictor.web_predictor([smi]) #has to be a list
		return res


#for testing purposes
if __name__ == "__main__":
	predictor = Site_Predictor()
	react = 'Cc1ccccc1'
	print(react)
	res = predictor.predict(react)
	print(res[0])

