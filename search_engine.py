class SearchEngineMeta(type):
    """Search engine metaclass that will be used for search engine class creation."""

    def __instancecheck__(cls, instance):
        return cls.__subclasscheck__(type(instance))

    def __subclasscheck__(cls, subclass):
        return hasattr(subclass, 'search') and \
               callable(subclass.search)


class SearchEngineInterface(metaclass=SearchEngineMeta):
    """This interface is used for concrete classes to inherit from."""
    pass
