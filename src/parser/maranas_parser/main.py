import csv

class Substrate:
    def __init__(self, size, id, name):
        self.size = size
        self.id = id
        self.name = name
        self.coefficient = 1.0

    def __str__(self):
        print(str(self.id) + ' ' + self.name + ' ', end='')
        return ""

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
        print('id: ' + str(self.id))
        print('name: ' + self.name)
        print('is_reversed: ' + str(self.is_reversed))
        print('is_excluded: ' + str(self.is_excluded))
        for sub in self.chemical_reaction.left:
            print(sub)
        for sub in self.chemical_reaction.right:
            print(sub)
        for transition in self.chemical_reaction.atom_transitions:
            print(transition.substrate_pos, transition.substrate_atom, transition.product_pos, transition.product_atom, sep=' ')
        print()
        return ""

def parse_chemical_equation_side(side, prefix=None):
    result = []
    is_excluded = False
    for substance in side.split(' '):
        if substance != '+' and substance:
            new_substance = Substrate(0, len(result), substance)
            if substance.startswith('0*'): # excluded metabolite
                new_substance.name = substance[2:]
                is_excluded = True
            if prefix:
                new_substance.name += prefix
            result.append(new_substance)
    return result, is_excluded


def parse_reaction(reaction, name):
    result = Reaction(name)

    separator = reaction.find('-->')
    separator_type = '-->'
    if separator == -1:
        separator = reaction.find('->')
        separator_type = '->'
        if separator == -1:
            separator = reaction.find('<==>')
            separator_type = '<==>'
            assert(separator != -1)

    left = reaction[0:separator]
    right = reaction[separator + len(separator_type):]
    prefix = None
    if left[0] == '[':
        prefix = left.split(' ')[0]
        prefix_end = left.find(':')
        assert(prefix_end != -1)
        left = left[prefix_end + 1:]

    left_side, is_excluded = parse_chemical_equation_side(left, prefix)
    right_side, _ = parse_chemical_equation_side(right, prefix)
    is_reversed = (separator_type == '<==>')
    result.is_reversed = is_reversed
    result.is_excluded = is_excluded
    result.chemical_reaction = ChemicalEquation(left_side, right_side)
    return result


if __name__ == '__main__':
    reactions = dict()
    with open('../../../modelMaranas/model.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        next(reader, None)
        for row in reader:
            atoms = row[0].split(',')
            substrate = row[2]
            side = row[3]
            reaction = row[4]
            name = row[5]
            if name not in reactions:
                reactions[name] = parse_reaction(reaction, name)

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
    for reaction in reactions.values():
        print(reaction)


