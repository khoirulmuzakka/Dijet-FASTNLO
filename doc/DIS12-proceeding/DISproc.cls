\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{DISproc}[2012/03/06 V1.1 DIS 2012 proceedings]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[twoside,a4paper]{article}
% \RequirePackage[dvips]{graphicx,epsfig,color}
\RequirePackage{graphicx}
\RequirePackage{wrapfig,rotating}
\RequirePackage{amssymb,amsmath,array}
\RequirePackage{booktabs}
\RequirePackage{hepunits}
\RequirePackage{csquotes}
% \RequirePackage[citecolor=magenta,linkcolor=blue,urlcolor=black]{hyperref}
\RequirePackage[urlcolor=black]{hyperref}

\setlength{\textwidth}{145mm}
\setlength{\textheight}{198mm}
\setlength{\topmargin}{5mm}
\setlength{\oddsidemargin}{6mm}
\setlength{\evensidemargin}{8mm}

\renewcommand{\heavyrulewidth}{0.1em}

\newif\if@doi\@doifalse
\def\contribID#1{\gdef\@contribID{#1}}
\def\@contribID{\@latex@warning{Please provide the Contribution-ID}}
\def\confID#1{\@doitrue\gdef\@confID{#1}}
\def\doi{\@doitrue}

\def\DISproc#1{\gdef\@DISproc{#1}}                       
\gdef\@DISproc{\@latex@warning{No DIS-PROC number given}}

\def\acronym#1{\gdef\@acronym{#1}}                                  
\let\@acronym\@empty                                     

% Conference information
% TO THE CONFERENCE EDITORS: 
% please update the following information      
% before sending the template to the authors
\confID{153252}  % if the conference is on Indico uncomment this line
\DISproc{DIS 2012 proceedings}
\acronym{DIS 2012} % if you want the acronym in the page footer uncomment this line


\def\@maketitle{%
  \newpage
  \parindent 0mm
  \null
  \vskip 2em%
  \let \footnote \thanks
    {\bfseries\LARGE \@title \par}%
    \vskip 1.5em%
    {\lineskip .5em%
        \@author
    }%
  \if@doi
     \par
     \vskip 1em
% =======  To the editor:  ==========================================
%   Please swap the comment of the following 2 lines. 
%   So the propper DOI appears on the document.
% -------------------------------------------------------------------
%     {\bfseries DOI:} {\texttt http://dx.doi.org/10.3204/\@DISproc/\@contribID}
     {\bfseries DOI:} {\texttt will be assigned}
  \fi
  \par
  \vskip 1.5em}

\if@titlepage
\else
  \renewenvironment{abstract}{%
      \if@twocolumn
        \section*{\abstractname}%
      \else
        \small
        \quote
      \fi}
      {\if@twocolumn\else\endquote\fi}
\fi

\def\ps@plain{\let\@mkboth\@gobbletwo
     \let\@oddhead\@empty\let\@evenhead\@empty
     \def\@evenfoot{\thepage\hfil\slshape \@acronym}%
     \def\@oddfoot{{\slshape \@acronym}\hfil\thepage}%
     }
\def\ps@empty{\let\@mkboth\@gobbletwo
     \let\@oddhead\@empty\let\@evenhead\@empty
     \def\@evenfoot{\hfil\slshape \@acronym}%
     \def\@oddfoot{{\slshape \@acronym}\hfil}%
}
\pagestyle{plain}
