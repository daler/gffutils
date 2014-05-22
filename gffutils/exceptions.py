class FeatureNotFoundError(Exception):
    """
    Error to be raised when an ID is not in the database.
    """
    def __init__(self, feature_id):
        Exception.__init__(self)
        self.feature_id = feature_id

    def __str__(self):
        return self.feature_id


class DuplicateIDError(Exception):
    pass


class AttributeStringError(Exception):
    pass
