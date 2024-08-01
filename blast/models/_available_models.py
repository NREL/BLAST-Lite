import sys, inspect


def available_models():
    """
    Return list of names of available models.
    """
    clsmembers = inspect.getmembers(sys.modules['blast.models'], inspect.isclass)
    return [cls[0] for cls in clsmembers if cls[0] != 'BatteryDegradationModel']