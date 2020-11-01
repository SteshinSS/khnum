import csv

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
        return str(self.id) + ' ' + str(self.coefficient) + ' ' + self.name + ' '

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
        self.atom_transitions = []

class Reaction:
    count = 0

    def __init__(self, name):
        self.id = Reaction.count
        Reaction.count += 1
        self.name = name
        self.chemical_reaction = None
        self.is_reversed = None
        self.is_excluded = None

    def __str__(self):
        res = ""
        res += 'id: ' + str(self.id) + '\n'
        res += 'name: ' + self.name + '\n'
        res += 'is_reversed: ' + str(self.is_reversed) + '\n'
        res += 'is_excluded: ' + str(self.is_excluded) + '\n'
        for sub in self.chemical_reaction.left:
            res += str(sub) + '\n'
        for sub in self.chemical_reaction.right:
            res += str(sub) + '\n'
        for transition in self.chemical_reaction.atom_transitions:
            res += str(transition.substrate_pos) + ' '
            res += str(transition.substrate_atom) + ' '
            res += str(transition.product_pos) + ' '
            res += str(transition.product_atom) + '\n'
        return res

def parse_chemical_equation_side(side, prefix=None, raw=None):
    result = []
    is_excluded = False
    last_coefficient = None
    for substance in side.split(' '):
        if substance != '+' and substance:
            if is_number(substance):
                assert(last_coefficient is None)
                last_coefficient = float(substance)
            else:
                substance_name = substance
                if substance.startswith('0*'): # excluded metabolite
                    substance_name = substance_name[2:]
                    is_excluded = True
                if prefix:
                    substance_name += prefix
                atoms = [item for item in raw.atoms if item[0] == substance_name]
                new_substance = Substrate(len(result), substance_name)
                if atoms:
                    if new_substance.size == None:
                        new_substance.size = len(atoms[0][1])
                    else:
                        assert new_substance.size == len(atoms[0][1])
                else:
                    new_substance.size = 0



                if last_coefficient:
                    new_substance.coefficient = last_coefficient
                if atoms:
                    new_substance.coefficient = 1.0
                    sub = new_substance
                    for i in range(len(atoms)):
                        result.append(sub)
                        sub = Substrate(len(result), substance_name)
                        sub.size = len(atoms[0][1])
                else:
                    result.append(new_substance)
                last_coefficient = None
    return result, is_excluded


def parse_reaction(raw):
    result = Reaction(raw.name)
    reaction1 = raw.reaction

    separator = reaction1.find('-->')
    separator_type = '-->'
    if separator == -1:
        separator = reaction1.find('->')
        separator_type = '->'
        if separator == -1:
            separator = reaction1.find('<==>')
            separator_type = '<==>'
            assert(separator != -1)

    left = reaction1[0:separator]
    right = reaction1[separator + len(separator_type):]
    prefix = None
    if left[0] == '[':
        prefix = left.split(' ')[0]
        prefix_end = left.find(':')
        assert(prefix_end != -1)
        left = left[prefix_end + 1:]

    left_side, is_excluded = parse_chemical_equation_side(left, prefix, raw)
    right_side, _ = parse_chemical_equation_side(right, prefix, raw)
    is_reversed = (separator_type == '<==>')
    result.is_reversed = is_reversed
    result.is_excluded = is_excluded
    result.chemical_reaction = ChemicalEquation(left_side, right_side)
    atom_transitions = []
    left_pos_to_fill = dict()
    right_pos_to_fill = dict()
    while True:
        if not raw.atoms:
            break
        sub, atoms = raw.atoms[0]
        if not atoms:
            raw.atoms.pop(0)
            continue

        is_left = False
        sub_position = None
        for sub_pos in range(len(left_side)):
            if left_side[sub_pos].name == sub:
                if sub_pos in left_pos_to_fill:
                    if left_pos_to_fill[sub_pos] < left_side[sub_pos].size:
                        is_left = True
                        left_pos_to_fill[sub_pos]+= 1
                        sub_position = sub_pos
                        break
                else:
                    left_pos_to_fill[sub_pos] = 1
                    is_left = True
                    sub_position = sub_pos
                    break

        if not is_left:
            for sub_pos in range(len(right_side)):
                if right_side[sub_pos].name == sub:
                    if sub_pos in right_pos_to_fill:
                        if right_pos_to_fill[sub_pos] < right_side[sub_pos].size:
                            right_pos_to_fill[sub_pos] += 1
                            sub_position = sub_pos
                            break
                    else:
                        right_pos_to_fill[sub_pos] = 1
                        sub_position = sub_pos
                        break


        atom = atoms.pop(0)
        this_atom = [pair for pair in raw.atoms if atom in pair[1]]
        assert len(this_atom) == 1
        prod = this_atom[0][0]
        this_atom[0][1].pop(this_atom[0][1].index(atom))

        prod_position = None

        if is_left:
            for prod_pos in range(len(right_side)):
                if right_side[prod_pos].name == prod:
                    if prod_pos in right_pos_to_fill:
                        if right_pos_to_fill[prod_pos] < right_side[prod_pos].size:
                            right_pos_to_fill[prod_pos] += 1
                            prod_position = prod_pos
                            break
                    else:
                        right_pos_to_fill[prod_pos] = 1
                        prod_position = prod_pos
                        break
        else:
            for prod_pos in range(len(left_side)):
                if left_side[prod_pos].name == prod:
                    if prod_pos in left_pos_to_fill:
                        if left_pos_to_fill[prod_pos] < left_side[prod_pos].size:
                            left_pos_to_fill[prod_pos] += 1
                            prod_position = prod_pos
                            break
                    else:
                        left_pos_to_fill[prod_pos] = 1
                        prod_position = prod_pos
                        break
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
    result.chemical_reaction.atom_transitions = atom_transitions
    return result


class RawReaction():
    def __init__(self):
        self.name = None
        self.reaction = None
        self.atoms = []

if __name__ == '__main__':
    raw_reactions = dict()
    reactions = dict()
    with open('../../../modelMaranas/model.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            atoms = row[0].split(',')
            substrate = row[2]
            reaction = row[4]
            name = row[5]
            if name not in raw_reactions:
                raw = RawReaction()
                raw.name = name
                raw.reaction = reaction
                raw_reactions[name] = raw

            new_atoms = [substrate, atoms]
            raw_reactions[name].atoms.append(new_atoms)

        for raw in raw_reactions.values():
            reaction = parse_reaction(raw)
            print(reaction)

            """
            parse_reaction(reaction, name)

            substrate_pos = None
            sub_in_left = [i for i, e in enumerate(reactions[name].chemical_reaction.left) if e.name == substrate]
            sub_in_right = [i for i, e in enumerate(reactions[name].chemical_reaction.right) if e.name == substrate]
            if sub_in_left:
                side = 'reactant'
            else:
                side = 'product'
            for atom in atoms:
                if side == 'reactant':
                    substrate_pos = [i for i, e in enumerate(reactions[name].chemical_reaction.left) if e.name == substrate][0]
                    transition = [i for i in reactions[name].chemical_reaction.atom_transitions if i.product_atom == atom]
                    if not transition:
                        transition = AtomTransition()
                        transition.substrate_pos = substrate_pos
                        transition.substrate_atom = atom
                        reactions[name].chemical_reaction.atom_transitions.append(transition)
                    else:
                        transition[0].substrate_pos = substrate_pos
                        transition[0].substrate_atom = atom
                elif side == 'product':
                    substrate_pos = [i for i, e in enumerate(reactions[name].chemical_reaction.right) if e.name == substrate][0]
                    transition = [i for i in reactions[name].chemical_reaction.atom_transitions if
                                  i.substrate_atom == atom]
                    if not transition:
                        transition = AtomTransition()
                        transition.product_pos = substrate_pos
                        transition.product_atom = atom
                        reactions[name].chemical_reaction.atom_transitions.append(transition)
                    else:
                        transition[0].product_pos = substrate_pos
                        transition[0].product_atom = atom
            """
    for reaction in reactions.values():
        print(reaction)


