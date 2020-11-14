def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class Substrate:
    def __init__(self, id, name):
        self.size = None
        self.id = id
        self.name = name
        self.coefficient = 1.0

    def __str__(self):
        return str(self.id) + ' ' + str(self.coefficient) + ' ' + self.name + ' ' + str(self.size) + ' '


class AtomTransition:
    def __init__(self):
        self.substrate_pos = None
        self.product_pos = None
        self.substrate_atom = None
        self.product_atom = None


class ChemicalEquation:
    def __init__(self, left, right):
        self.left = left
        self.right = right
        self.atom_transitions = []  # contains AtomTransitions


class Reaction:
    count = 0

    def __init__(self, name):
        self.is_flux = False
        self.id = Reaction.count
        Reaction.count += 1
        self.name = name
        self.chemical_reaction = None
        self.is_reversed = None
        self.is_excluded = None  # for reactions with 0* substrates

    def __str__(self):  # is needed for string forming
        res = ""
        res += str(self.id) + '\n'
        res += self.name + '\n'

        res += str(int(self.is_flux)) + '\n'
        res += str(int(self.is_reversed)) + '\n'
        res += str(int(self.is_excluded)) + '\n'
        res += str(len(self.chemical_reaction.left)) + '\n'
        for sub in self.chemical_reaction.left:
            res += str(sub) + '\n'
        res += str(len(self.chemical_reaction.right)) + '\n'
        for sub in self.chemical_reaction.right:
            res += str(sub) + '\n'
        res += str(len(self.chemical_reaction.atom_transitions)) + '\n'
        for transition in self.chemical_reaction.atom_transitions:
            res += str(transition.substrate_pos) + ' '
            res += str(transition.substrate_atom) + ' '
            res += str(transition.product_pos) + ' '
            res += str(transition.product_atom) + '\n'
        return res


def is_sides_equals(first: [], second: []):
    import math
    if len(first) != len(second):
        return False

    is_found = [False] * len(second)
    for sub in first:
        here = False
        for i in range(len(second)):
            if not is_found[i]:
                if second[i].name == sub.name and math.isclose(sub.coefficient, second[i].coefficient):
                    is_found[i] = True
                    here = True
                    break
        if not here:
            return False
    return True

def is_reaсtions_equal(first: Reaction, second: Reaction):
    normal = is_sides_equals(first.chemical_reaction.left, second.chemical_reaction.left) and is_sides_equals(first.chemical_reaction.right, second.chemical_reaction.right)
    if normal:
        return True
    rev = is_sides_equals(first.chemical_reaction.left, second.chemical_reaction.right) and is_sides_equals(first.chemical_reaction.right, second.chemical_reaction.left)
    return rev

def parse_chemical_equation_side(side, prefix=None, raw_atoms=None):
    result = []
    is_excluded = False
    last_coefficient = None
    for substance in side.split(' '):
        if substance == '+' or not substance:
            continue

        if is_number(substance):
            assert(last_coefficient is None)  # two coefficients in a row
            last_coefficient = float(substance)
        else:
            substance_name = substance

            if '[' in substance_name:
                start = substance_name.index('[')
                stop = substance_name.index(']')
                pre = substance_name[start:stop + 1]

            if substance.startswith('0*'):  # excluded metabolite
                substance_name = substance_name[2:]
                is_excluded = True
            if prefix:
                substance_name += prefix
            atoms = [item for item in raw_atoms if item[0] == substance_name]
            new_substance = Substrate(len(result), substance_name)
            if last_coefficient:
                new_substance.coefficient = last_coefficient

            """ In case of several molecules of a substrate with a tracer's atoms are used in reaction, 
            we will write each molecule as separate. So instead of 4 ATP we form reaction: ATP + ATP + ATP + ATP."""

            if atoms:
                if new_substance.size is None:
                    new_substance.size = len(atoms[0][1])
                else:
                    assert new_substance.size == len(atoms[0][1])

                new_substance.coefficient = 1.0
                sub = new_substance
                for i in range(len(atoms)):
                    result.append(sub)
                    sub = Substrate(len(result), substance_name)
                    sub.size = len(atoms[0][1])
            else:
                new_substance.size = 0
                result.append(new_substance)
            last_coefficient = None
    return result, is_excluded


def get_reaction_separator(reaction):
    separator_pos = reaction.find('-->')
    separator_type = '-->'
    if separator_pos == -1:
        separator_pos = reaction.find('->')
        separator_type = '->'
        if separator_pos == -1:
            separator_pos = reaction.find('<==>')
            separator_type = '<==>'
            assert (separator_pos != -1)
    return separator_pos, separator_type


def find_substrate(side, substrate_name, pos_to_fill):
    is_found = False
    position = None
    for sub_pos in range(len(side)):
        if side[sub_pos].name == substrate_name:
            # if we found any atom transitions before
            if sub_pos in pos_to_fill:
                # it's possible to have several substrates with the same name in reaction
                # we will take first for which we don't know all atom transitions
                if pos_to_fill[sub_pos] < side[sub_pos].size:
                    is_found = True
                    pos_to_fill[sub_pos] += 1
                    position = sub_pos
                    break
            else:
                pos_to_fill[sub_pos] = 1
                is_found = True
                position = sub_pos
                break
    return is_found, position


def check_atom_transitions(side, pos_to_fill):
    for pos in range(len(side)):
        if side[pos].size > 0:
            assert pos_to_fill[pos] == side[pos].size


def get_atom_transitions(raw, left_side, right_side):
    atom_transitions = []

    # Maps substrate position to number of already processed atom transitions for it
    left_pos_to_fill = dict()
    right_pos_to_fill = dict()
    while True:
        if not raw.atoms:
            break
        substrate, atoms = raw.atoms[0]
        if not atoms:
            raw.atoms.pop(0)
            continue
        is_left, sub_position = find_substrate(left_side, substrate, left_pos_to_fill)
        if not is_left:
            _, sub_position = find_substrate(right_side, substrate, right_pos_to_fill)

        # looking for this atom transition
        atom = atoms.pop(0)
        this_atom = [pair for pair in raw.atoms if atom in pair[1]]
        assert len(this_atom) == 1
        product = this_atom[0][0]
        this_atom[0][1].pop(this_atom[0][1].index(atom))

        prod_position = None

        if is_left:
            is_found, prod_position = find_substrate(right_side, product, right_pos_to_fill)
        else:
            is_found, prod_position = find_substrate(left_side, product, left_pos_to_fill)
        assert is_found

        transition = AtomTransition()
        transition.substrate_pos = sub_position
        if is_left:
            transition.substrate_atom = left_pos_to_fill[sub_position] - 1
        else:
            transition.substrate_atom = right_pos_to_fill[sub_position] - 1
        transition.product_pos = prod_position
        if is_left:
            transition.product_atom = right_pos_to_fill[prod_position] - 1
        else:
            transition.product_atom = left_pos_to_fill[prod_position] - 1
        atom_transitions.append(transition)

    check_atom_transitions(left_side, left_pos_to_fill)
    check_atom_transitions(right_side, right_pos_to_fill)

    return atom_transitions


def parse_reaction(raw):
    result = Reaction(raw.name)
    reaction = raw.reaction
    separator_pos, separator_type = get_reaction_separator(reaction)

    left = reaction[0:separator_pos]
    right = reaction[separator_pos + len(separator_type):]
    prefix = None
    if left[0] == '[':
        prefix = left.split(' ')[0]
        prefix_end = left.find(':')
        assert(prefix_end != -1)
        left = left[prefix_end + 1:]

    left_side, is_excluded = parse_chemical_equation_side(left, prefix, raw.atoms)
    right_side, _ = parse_chemical_equation_side(right, prefix, raw.atoms)
    is_reversed = (separator_type == '<==>')
    result.is_reversed = is_reversed
    result.is_excluded = is_excluded
    result.chemical_reaction = ChemicalEquation(left_side, right_side)
    result.chemical_reaction.atom_transitions = get_atom_transitions(raw, left_side, right_side)
    return result


class RawReaction():
    def __init__(self):
        self.name = None
        self.reaction = None
        self.atoms = []  # contains pairs [substrate name, it's atoms]


def check_path_exists(path):
    import os
    if not os.path.exists(path):
        print("No file found on: " + path)
        return False
    else:
        return True


# Forms a string which will be parsed by C++
def print_reactions(reactions):
    result = str(len(reactions)) + '\n'
    for reaction in reactions.values():
        result += str(reaction) + '\n'
    return result.encode('utf-8')



""" 
The Parser itself. Gets path, return string for further parsing in such format:
id
name
is reaction reversed?
is there are excluded substrate? (eg 0*val-L[c] -> Val[d])
number of substrates in the left side of reaction
[id coefficient name number_of_atoms]
number of substrates in the right side of reaction
[id coefficient name number_of_atoms]
number of atom transitions
[substrate_position substrate_atom product_position product_atom]
"""
def parse(path):
    assert check_path_exists(path)

    emu_model = path + 'model.csv'

    raw_reactions = dict()
    reactions = dict()
    import csv
    csvfile = open(emu_model, 'r', newline='')
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(reader, None)  # skip header

    # first collect all atom transitions in one structure
    for row in reader:
        reaction = row[4]
        name = row[5]
        if name not in raw_reactions:
            raw = RawReaction()
            raw.name = name
            raw.reaction = reaction
            raw_reactions[name] = raw

        atoms = row[0].split(',')
        substrate = row[2]
        new_atoms = [substrate, atoms]
        raw_reactions[name].atoms.append(new_atoms)

    # now parse all atom transitions at once
    for raw in raw_reactions.values():
        reactions[raw.name] = parse_reaction(raw)

    flux_model = path + 'flux_model.csv'
    OnlyFlux = ["PPA", "DPCOAK", "NADK", "NADS1", "SULRi", "DM_4HBA", "DM_HMFURN", "biomass_out", "EX_ca2(e)",
                "EX_cl(e)", "EX_cobalt2(e)", "EX_cu2(e)", "EX_fe2(e)", "EX_h(e)", "EX_h2(e)", "EX_h2o(e)", "EX_k(e)",
                "EX_mg2(e)", "EX_mn2(e)", "EX_mobd(e)", "EX_nh4(e)", "o2_in", "EX_pi(e)", "EX_so4(e)", "EX_zn2(e)",
                "CAt6pp", "CLt3_2pp", "COBALT2tpp", "CU2tpp", "FE2tpp", "FEROpp", "Kt2pp", "MG2tpp", "MN2t3pp",
                "MN2tpp", "NAt3_1p5pp", "NAt3_2pp", "NAt3pp", "NH4tpp", "NI2t3pp", "NI2tpp", "O2tpp", "PIt2rpp",
                "ZN2t3pp", "ZN2tpp", "CYTBD2pp", "CYTBDpp", "CYTBO3_4pp", "NADH10", "NADH16pp", "NADH17pp", "NADH5",
                "NADPHQR2", "NADPHQR3", "NADTRHD", "THD2pp", "TRDR", "H2Otpp", "H2tpp", "MNt2pp", "H2Otex", "CA2tex",
                "CLtex", "COBALT2tex", "CU2tex", "FE2tex", "FE3tex", "H2tex", "Htex", "Ktex", "MG2tex", "MNtex",
                "MOBDtex", "NH4tex", "O2tex", "PItex", "SO4tex", "Zn2tex", "CAT", ]

    csvfile = open(flux_model, 'r', newline='')
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for row in reader:
        name = row[0]
        if name in OnlyFlux:
            reaction = row[2]
            raw = RawReaction()
            raw.name = name + '_flux'
            raw.reaction = reaction
            flux = parse_reaction(raw)
            flux.is_flux = True
            reactions[name + '_flux'] = flux
    return print_reactions(reactions)




if __name__=='__main__':
    result = parse('../../../modelMaranas/')
    # print(result.decode('utf-8'))


