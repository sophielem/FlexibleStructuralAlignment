#!/usr/bin/env python

"""
reindex_pdb.py startindex infile.pdb outfile.pdb
    Rearrange residue number of "infile.pdb" so the first residue start from
    "startindex". if the structure contains any missing residue, residue
    number gaps will be removed.

reindex_pdb.py seq.fasta  infile.pdb outfile.pdb
    Rearrange residue number of "infile.pdb" according to "seq.fasta" so the
    residue index is the same as that in "seq.fasta"

@author : zhanglab
"""


def reindex_pdb_by_index(startindex, pdb_txt):
    """
    Reindex residue number of PDB format text

    Args:
        startindex: index of first residue
        pdb_txt:    text of input PDB to be reindexed

    Returns:
        The pdb reindexed
    """
    pdb_txt_reindex = ''
    current_old_index = ''  # residue number in origin PDB
    warn_chain_id = ''  # warning about new chain ID
    res_seq_new = ''
    for line in pdb_txt.splitlines():
        if (len(line) < 27 or
                (not line.startswith("ATOM  ") and
                 not line.startswith("HETATM") and not line.startswith("TER"))):
            pdb_txt_reindex += line[:22] + res_seq_new + '\n'
            continue
        elif not line[16] in ['A', ' ']:  # alternative location identifier
            continue
        res_seq = line[22:27]  # residue sequence number
        current_chain_id = line[21]  # chain identifier

        if not current_old_index:  # first residue encountered
            current_old_index = res_seq  # residue number in origin PDB
            current_new_index = int(startindex)
            chain_id = current_chain_id
            res_seq_new = str(current_new_index)
            res_seq_new = ' '*(4-len(res_seq_new)) + res_seq_new + ' '
        elif current_chain_id != chain_id:
            if warn_chain_id != current_chain_id and not line.startswith("TER"):
                return None
        elif res_seq != current_old_index:
            current_new_index += 1
            current_old_index = res_seq
            res_seq_new = str(current_new_index)
            res_seq_new = ' '*(4-len(res_seq_new)) + res_seq_new + ' '
        pdb_txt_reindex += (line[:16] + ' ' + line[17:22] +
                            res_seq_new + line[27:] + '\n')
    return pdb_txt_reindex


def reindex_pdb(startindex, infile):
    """
    parse PDB file "infile", reindex it according to start index or
    sequence file "startindex", and return the text of renumbered PDB

    Args:
        infile: Input pdb file to parse
        startindex: Start index to reindex

    Returns:
        The pdb reindexed
    """
    clean = True
    filin = open(infile, 'r')
    pdb_txt = ''
    for line in filin.read().splitlines():
        if line.startswith("END"):
            line = line.replace("ENDMDL", "END   ")
            pdb_txt += line+'\n'
            break
        if (line.startswith("ATOM  ") or line.startswith("TER") or
                (not clean and not line[:6] in ["DBREF ", "SEQADV", "MODRES",
                                                "HELIX ", "SHEET ", "SSBOND",
                                                "SITE  "])):
            pdb_txt += line + '\n'
    filin.close()

    pdb_txt_reindex = reindex_pdb_by_index(startindex, pdb_txt)
    return pdb_txt_reindex


def re_number_pdb(start, input_pdb, output_pdb):
    """
    parse PDB file, reindex it according to start index given,
    and write the new pdb reindexed

    Args:
        start: Start index to reindex
        input_pdb: Input pdb file to parse
        output_pdb: Output pdb reindexed

    Returns:
        True if all is correctly done or False if the pdb has several chains
    """
    # parse PDB file
    pdb_txt_reindex = reindex_pdb(start, input_pdb)

    # write PDB file
    if pdb_txt_reindex is not None:
        filin = open(output_pdb, 'w')
        filin.write(pdb_txt_reindex)
        filin.close()
        return True

    return False
