import re
import os

def main():
    """
    Automatically generates a beamer presentation draft from an article
    Notes:
           - Some functionality has been disabled by comments.
    """

    article_file_path = 'mek4420-eksamen-Ole.tex'
    beamer_file_path = 'slides.tex'

    written_figures = []

    with open(article_file_path, 'r') as art:
        with open(beamer_file_path, 'w') as outf:
           # in_frame = False
            for line in art:

                # Preamble
                #(currently ignoring)

                # Copy header
                # Ignores unnumbered sections
                if line.startswith(r'\section{') or line.startswith(r'\subsection'):

                 #   if in_frame:
                  #      outf.write(r'\end{frame}'+'\n')
                   #     in_frame = False

                    outf.write(safe_line(line))

                if line.startswith(r'\subsubsection'):
                    outf.write('%' + safe_line(line))

                    # Create frame for subsections
               #     if line.startswith(r'\subsection'):
                #        outf.write(r'\begin{frame}{\secname: \subsecname}'+'\n')
                 #       in_frame = True

                # Insert figure
                fig_refs = re.findall(r'\\ref{fig:([\w\d -_]+)}', line)
                if fig_refs:
                    for fig_ref in fig_refs:
                        if fig_ref not in written_figures:
                            fig_frame = safe_line(find_fig_str(article_file_path, fig_ref))
                            outf.write(fig_frame)
                            written_figures.append(fig_ref)

            #if in_frame:
             #   outf.write(r'\end{frame}' + '\n')
              #  in_frame = False

def safe_line(line):
    """ Processes latex control sequences which cause beamer to break
     * Remove \label{}, \todo{}
     * Add \protect infront of \textit{}
    """
    sline = line
    sline = re.sub(r'\\label{[\w\d:\s]+}', '', sline)
    sline = re.sub(r'\\todo{[\w\d:\s]+}', '', sline)

    sline = re.sub(r'\\textit{', r'\\protect\\textit{', sline)
    return sline

def find_fig_str(text_file_path, fig_ref):
    """ Finds the figure floating environment based on label
    """
    with open(text_file_path, 'r') as tfile:
        envstr = ''
        in_fig = False
        for line in tfile:

            if line.startswith(r"\begin{figure}"):
                in_fig = True

            if in_fig:
                to_write = line
                # Remove placement specifications:
                if line.startswith(r"\begin{figure}"):
                    to_write = re.sub(r'\[([\w]+)\]', '', line)
                # ignore label
                if line.startswith(r'\label{'):
                    to_write = ''
                # replace figure width to make it fill the whole slide
                if re.match(r'[\t\s]?\\includegraphics',line):
                    to_write = re.sub(r'\[([\w=\d]+)\]', r'[width =\\textwidth, height = 0.6\\textheight, keepaspectratio]', line)

                envstr += to_write

            labelmatch = re.match(r'[\t\s]?\\label{fig:([\w\d -_]+)}',line)
            if labelmatch:
                if labelmatch.group(1) == fig_ref:
                    right_fig = True
                else:
                    right_fig = False

            if line.startswith(r"\end{figure}"):
                in_fig = False
                if right_fig: # Stop searching
                    break
                else:
                    envstr = ''

    # Make figure be in its own frame:
    if envstr:
        envstr = r'\begin{frame}{\subsecname}' + '\n' + envstr + r'\end{frame}' + '\n'

    return envstr

main()