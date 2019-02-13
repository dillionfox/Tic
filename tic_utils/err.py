
class Error(Exception):
	"""
	Base class to handle custom exceptions

	"""
	pass

class InputError(Error):
	"""
	Custom error messages for mis-handling inputs

	"""
	def __init__(self, expr, message):
		self.expr = expr
		self.message = message

