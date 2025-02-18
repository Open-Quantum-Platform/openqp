"""Input parser"""
import configparser

class OQPConfigParser(configparser.ConfigParser):
    """Extends configparser with validation schema"""
    def __init__(self, *args, schema=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.schema = schema
        if self.schema:
            self.read_dict(self.strip_schema())

    def load_dict(self, input_dict):
        """Assign dict to config"""
        for section in input_dict.keys():
            for option, value in input_dict[section].items():
                self[section][option] = value

    def print_config(self):
        """Print resulting config"""
        for section in self.sections():
            for option, value in self[section].items():
                print(f'{section}.{option}={value}')

    def strip_schema(self):
        """Convert input schema to dictionary used by OQP"""
        stripped = {}
        for section, options in self.schema.items():
            stripped[section] = {}
            for name, descr in options.items():
                value = descr['default']
                stripped[section][name] = value
        return stripped

    def validate(self):
        """Validate configuration"""
        config = {}
        if not self.has_section('input'):
            raise ValueError("Missing section [input]")

        if not self.has_option('input', 'system'):
            raise ValueError("Missing molecule description")

        for section in self.sections():
            if section not in self.schema:
                raise ValueError(f"Unknown section: {section}")

            config[section] = {}
            if section == 'tests':
                for option, value in self[section].items():
                    config[section][option] = value
            else:
                valid_options = self.schema[section].keys()
                for option, value in self[section].items():
                    if option not in valid_options:
                        raise ValueError(f"Unknown option: {section}.{option}")
                    converter = self.schema[section][option]['type']
                    if converter == bool:
                        config[section][option] = self.getboolean(section, option)
                    else:
                        config[section][option] = converter(value)

        return config
