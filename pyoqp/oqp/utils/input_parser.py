"""Input parser"""
import configparser
import io

class OQPConfigParser(configparser.ConfigParser):
    """Extends configparser with validation schema"""
    def __init__(self, *args, schema=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.optionxform = str
        self.schema = schema
        if self.schema:
            self.read_dict(self.strip_schema())

    def _preprocess_input_file(self, filename):
        """Preprocess input file to handle system coordinates properly"""
        with open(filename, 'r') as f:
            lines = f.readlines()

        processed_lines = []
        in_system = False
        current_system = None
        system_coords = []

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith('['):
                if in_system and system_coords:
                    processed_lines.append('\n'.join(system_coords))
                in_system = False
                current_system = None
                processed_lines.append(line)
                continue

            if line.startswith('system'):
                if '=' in line:  
                    if in_system and system_coords:
                        processed_lines.append('\n'.join(system_coords))
                        system_coords = []
                    in_system = True
                    current_system = line.split('=')[0].strip()
                    processed_lines.append(line)
                    continue

            if in_system:
                if any(line.startswith(k + '=') for k in ['charge', 'runtype', 'basis', 'functional', 'method', 'd4']):
                    in_system = False
                    if system_coords:
                        processed_lines.append('\n'.join(system_coords))
                        system_coords = []
                    current_system = None
                    processed_lines.append(line)
                else:
                    if not line.startswith(' '):
                        line = '        ' + line
                    system_coords.append(line)
            else:
                processed_lines.append(line)

        if in_system and system_coords:
            processed_lines.append('\n'.join(system_coords))

        return '\n'.join(processed_lines)

    def read(self, filenames, encoding=None):
        """Override read to preprocess input files"""
        if isinstance(filenames, str):
            filenames = [filenames]
        
        for filename in filenames:
            processed_content = self._preprocess_input_file(filename)
            super().read_string(processed_content)

        return self

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
