import subprocess

class OPSIN_Lookup:
	def __init__(self):
		self.opsin_path = 'opsin-2.0.0-jar-with-dependencies.jar'
		self.proc = subprocess.Popen(['java', '-jar', self.opsin_path], shell = False, stdin = subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		self.connected = True

	def lookup(self, name):
		self.proc.stdin.write('%s\n' % name)
		return self.proc.stdout.readline().strip()

	def __del__(self):
		self.proc.stdin.close()
		self.proc.wait()
		self.connected = False