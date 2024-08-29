"""
Molden compatibility module
reader class
"""

import re
from oqp.utils.matrix import *

class MoldenReader:
    """
     Class to read and process data from Molden format
     Reads alpha and beta MOs vectors and energies
    """

    def __init__(self, filename):
        self.sections = None
        self.filename = filename

    def _extract_sections_from_file(self):
        """Extracts sections of data from Molden file"""
        with open(self.filename, 'r') as file:
            content = file.read()

        section_pattern = r'\[(\w+)\]\s*((?:(?!\[\w+\]).)*)'
        matches = re.findall(section_pattern, content, re.DOTALL)

        sections_dict = {}
        for section_name, section_content in matches:
            sections_dict[section_name] = section_content.strip()

        if not sections_dict["MO"]:  # Check if the "MO" section is empty
            raise ValueError("[MO] section is empty in the file.")

        return sections_dict

    @staticmethod
    def _parse_molecular_orbitals(mo_data):
        """Parses Molecular Orbitals data"""
        coef_alpha, coef_beta, energy_alpha, energy_beta = [], [], [], []
        occupancy_alpha, occupancy_beta = [], []

        for entry in mo_data:
            if entry["spin"] == 'Alpha':
                coef_alpha.append(entry["coefficients"])
                energy_alpha.append(entry["energy"])
                if entry["occupancy"]:
                    occupancy_alpha.append(entry["occupancy"])
            elif entry["spin"] == 'Beta':
                coef_beta.append(entry["coefficients"])
                energy_beta.append(entry["energy"])
                if entry["occupancy"]:
                    occupancy_beta.append(entry["occupancy"])

        density_alpha = orb_to_dens(np.array(coef_alpha), np.array(occupancy_alpha))
        res = {
            "mo_vec_a": coef_alpha,
            "mo_e_a": energy_alpha,
            "dens_a": density_alpha,
        }

        if any(entry["spin"] == 'Beta' for entry in mo_data):
            density_beta = orb_to_dens(np.array(coef_beta), np.array(occupancy_beta))
            res["mo_vec_b"] = coef_beta
            res["mo_e_b"] = energy_beta
            res["dens_b"] = density_beta

        return res

    def read_mo(self):
        """Extracts and returns Molecular Orbital data from the file"""

        try:
            # Try to access the "MO" section
            self.sections = self._extract_sections_from_file()
        except KeyError:
            # KeyError will be raised if the "MO" key is not found in self.sections
            raise KeyError("[MO] section is missing or malformed in the file.")

        individual_mo_data = re.split(r'\n(?=Ene=)', self.sections["MO"])

        parsed_data = []
        for mo_section in individual_mo_data:
            lines = mo_section.strip().splitlines()
            mo_details = {}
            for line in lines:
                if line.startswith("Ene="):
                    mo_details["energy"] = float(line.split("=")[1].strip())
                elif line.startswith("Spin="):
                    mo_details["spin"] = line.split("=")[1].strip()
                elif line.startswith("Occup="):
                    mo_details["occupancy"] = float(line.split("=")[1].strip())
                else:
                    try:
                        _, coef = line.split()
                        mo_details.setdefault("coefficients", []).append(float(coef))
                    except ValueError:
                        pass

            required_keys = {"energy", "spin", "occupancy", "coefficients"}
            if required_keys.issubset(mo_details.keys()):
                parsed_data.append(mo_details)

        return self._parse_molecular_orbitals(parsed_data)
