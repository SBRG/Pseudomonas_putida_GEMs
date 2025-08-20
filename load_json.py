import json
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, Union
from cobra import Gene, Metabolite, Model, Reaction
from cobra.util import set_objective
from collections import OrderedDict

def _metabolite_from_dict(metabolite):
    """Convert a dictionary to cobra Metabolite object.

    Parameters
    ----------
    metabolite : dict
        The dictionary to convert to cobra.Metabolite .

    Returns
    -------
    cobra.Metabolite

    See Also
    --------
    _metabolite_to_dict : Convert a cobra Metabolite object to dictionary.

    """
    new_metabolite = Metabolite()
    for k, v in metabolite.items():
        setattr(new_metabolite, k, v)
    return new_metabolite

def _reaction_from_dict(reaction, model):
    """Convert a dictionary to a cobra Reaction object.

    Parameters
    ----------
    reaction : dict
        The dictionary to convert to cobra.Reaction .
    model : cobra.Model
        The model to which the reaction should associate with.

    Returns
    -------
    cobra.Reaction
        The converted cobra.Reaction object.

    See Also
    --------
    _reaction_to_dict : Convert a cobra Reaction object to a dictionary.

    """
    new_reaction = Reaction()
    for k, v in reaction.items():
        if k in {"objective_coefficient", "reversibility", "reaction"}:
            continue
        elif k == "metabolites":
            new_reaction.add_metabolites(
                OrderedDict(
                    (model.metabolites.get_by_id(str(met)), coeff)
                    for met, coeff in v.items()
                )
            )
        else:
            if k == "lower_bound" or k == "upper_bound":
                setattr(new_reaction, k, float(v))
            else:
                setattr(new_reaction, k, v)
    return new_reaction


def gene_from_dict(gene):
    """Convert a dictionary to cobra Gene object.

    Parameters
    ----------
    gene : dict
        The dictionary to convert to cobra.Gene .

    Returns
    -------
    cobra.Gene
        The converted cobra.Gene object.

    See Also
    --------
    _gene_to_dict : Convert a cobra Gene object to a dictionary.

    """
    new_gene = Gene(gene["id"])
    for k, v in gene.items():
        setattr(new_gene, k, v)
    return new_gene



def model_from_dict(obj):
    """Build a cobra Model from a dictionary.

    Models stored in JSON are first formulated as a dictionary that can be read
    to a cobra Model using this function.

    Parameters
    ----------
    obj : dict
        A dictionary with keys: 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists of dictionaries holding all
        attributes to form the corresponding object.

    Returns
    -------
    cobra.Model
        The generated model.

    Raises
    ------
    ValueError
        If `obj` has no 'reactions' attribute.

    See Also
    --------
    model_to_dict : Convert a cobra Model to a dictionary.

    """
    if "reactions" not in obj:
        raise ValueError("Object has no .reactions attribute. Cannot load.")
    model = Model()
    model.add_metabolites(
        [_metabolite_from_dict(metabolite) for metabolite in obj["metabolites"]]
    )
    model.genes.extend([gene_from_dict(gene) for gene in obj["genes"]])
    model.add_reactions(
        [_reaction_from_dict(reaction, model) for reaction in obj["reactions"]]
    )
    objective_reactions = [
        rxn for rxn in obj["reactions"] if rxn.get("objective_coefficient", 0) != 0
    ]
    coefficients = {
        model.reactions.get_by_id(rxn["id"]): rxn["objective_coefficient"]
        for rxn in objective_reactions
    }
    set_objective(model, coefficients)
    for k, v in obj.items():
        if k in {"id", "name", "notes", "compartments", "annotation","subsystem"}:
            setattr(model, k, v)
    return model




def load_json_model(filename):
    """Load a cobra model from a file in JSON format.

    Parameters
    ----------
    filename : str or file-like
        File path or descriptor that contains the JSON document describing the
        cobra model.

    Returns
    -------
    cobra.Model
        The cobra model as represented in the JSON document.

    See Also
    --------
    from_json : Load from a JSON string.

    """
    if isinstance(filename, (str, Path)):
        with open(filename, "r") as file_handle:
            return model_from_dict(json.load(file_handle))
    else:
        return model_from_dict(json.load(filename))