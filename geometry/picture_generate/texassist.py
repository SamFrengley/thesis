import os
import pathlib


def format_Fg(*labs, lam=""):
    """
    Formats a string for the modular curve F_g^+.

    Parameters
    ----------
    *labs : {"I", "c", "ns", "s", "w"}
        Parts of the subscript for g.
    lam : str, default=""
        Optional multiply the parentheses by a constant.

    Returns
    -------
    str
        A string which nicely formats as LaTex code the curve F_g^+.
    """
    new = []
    for l in labs:
        if l == "I" or len(l) > 2:
            nl = l
        else:
            nl = f"\\mathrm{{{l}}}"
        new.append(nl)

    if len(labs) == 1:
        ret = f"$\\widetilde{{F}}_{{{lam}{new[0]}}}^+$"
    else:
        this_lab = ",".join(new)
        ret = f"$\\widetilde{{F}}_{{{lam}({this_lab})}}^+$"
    return ret


def careful_dir_name(rel_dir_name):
    """
    Strips any trailing forward slashes from a string.

    Parameters
    ----------
    rel_dir_name : str
        A string which will look like *foo/bar/ or *foo/bar.

    Returns
    -------
    str
        `rel_dir_name` but with trailing / removed
    """
    if rel_dir_name == "":
        rel_dir_name = "."
    elif rel_dir_name[-1] == "/":
        rel_dir_name = rel_dir_name[:-1]
    return rel_dir_name


def make_tex_file(rel_dir_name, pgf_filenames, captions=[]):
    """
    Makes two tex files in a relative directory which import all listed pgf 
    files and gives optional captions to the corresponding figures. The files 
    are:
        main.tex : When run produces a pdf with the images.
        source-pgfs.tex : The functional part of main.tex, can be imported into 
            other .tex files

    Parameters
    ----------
    rel_dir_name : str
        The directory in which the pgf files are, and in which to create the 
        *.tex files.
    pgf_filenames : list
        A list of pgf files which are to be imported in a .tex file
    captions : list, default=[]
        A list captions which if non-empty should have the same length as 
        pfg_filenames and gives the figures the captions listed.

    Returns
    -------
    None
    """
    main_str = "\\documentclass{article}\n"
    main_str += "\\usepackage[a4paper, margin=2cm]{geometry}\n"
    main_str += "\\usepackage{pgf}\n"
    main_str += "\\pagestyle{empty}\n\n"
    main_str += "\\begin{document}\n"
    main_str += "\\input{source-pgfs.tex}\n"
    main_str += "\\end{document}"

    source_pgfs_str = ""
    for i in range(0, len(pgf_filenames)):
        if len(captions) == 0:
            cpt = ""
        elif len(captions) == len(pgf_filenames):
            cpt = f"  \\caption{{{captions[i]}}}\n"
        else:
            raise ValueError("captions and pgf_filenames different lengths")
        source_pgfs_str += "\\begin{figure}[p]\n  \\centering\n"
        source_pgfs_str += f"  \\input{{{pgf_filenames[i]}}}\n"
        source_pgfs_str += cpt
        source_pgfs_str += "\\end{figure}\n\n"

    rel_dir_name = careful_dir_name(rel_dir_name)
    main_name = pathlib.Path(f"{rel_dir_name}/main.tex")
    source_pgfs_name = main_name.with_name("source-pgfs.tex")
    #    main_name = str(main_name)
    #    source_pgfs_name = str(source_pgfs_name)
    with main_name.open("w") as f:
        f.write(main_str)
    with source_pgfs_name.open("w") as f:
        f.write(source_pgfs_str)


def compile_tex_file(dir_name, filename="main.tex", open_pdf=False):
    """
    Compiles a tex file in a relative directory as if you were there.

    dir_name/main.tex as if you were there
    """
    dir_name = careful_dir_name(dir_name)
    filename = filename.strip()
    filename = pathlib.Path(filename)
    filename = filename.with_suffix(".tex")
    filename = str(filename)
    if not os.path.isfile(f"{dir_name}/{filename}"):
        raise FileNotFoundError(f"{filename} is not a file in {dir_name}")

    cd_str = f"cd {dir_name}"
    latex_str = f"pdflatex -interaction=batchmode {filename}"
    os.system(f"{cd_str} && {latex_str} && cd -")
    if open_pdf:
        filename = pathlib.Path(f"{dir_name}/{filename}")
        filename = filename.with_suffix(".pdf")
        os.system(f"open {str(filename)}")