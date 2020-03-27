
class Gap(object):
    def __init__(self, character="-"):
        self.character = "-"

    def __str__(self):
        return self.character

    def __eq__(self, other):
        return isinstance(other, Gap)