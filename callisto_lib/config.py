import os
import yaml
import traceback
import warnings

# the path to the default configuration file
_dirname = os.path.pardir(os.path.dirname(__file__))
_DEFAULT_CONFIG_FILE = os.path.join(_dirname, "default_config.yaml")

def load_config(filename = None):
    if filename is None:
        filename = _DEFAULT_CONFIG_FILE
    
    with open(filename, 'r') as stream:
        config = yaml.safe_load(stream)
        return config
    
    _config_file = os.environ.get("CALLISTO_CONFIG")
    
    if _config_file is None:
        config = load_config(_DEFAULT_CONFIG_FILE)
    else:
        try:
            config = load_config(_config_file)
        except Exception:
            msg = (
                "\nLoading configuration from CALLISTO_CONFIG enviromental "
                "variable failed:"
                "\n--- START IGNORED TRACEBACK --- \n"
                + traceback.format_exc()
                + "\n --- END IGNORED TRACEBACK ---"
                "\nLoading default configuration"
            )
            warnings.warn(msg)
            config = load_config(_DEFAULT_CONFIG_FILE)
    return