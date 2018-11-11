% -------------------------------------------------------------------
%  @LaTeX-class-file{
%     filename        = "DissertateUSU.cls",
%     version         = "1.0",
%     date            = "01 July 2017",
%     codetable       = "ISO/ASCII",
%     keywords        = "LaTeX, Dissertate",
%     note            = "Based on 'Dissertate' by Jordan Suchow,
%                        For support, email 'tyson.barrett@usu.edu',
%     docstring       = "Class for a dissertation at Utah State."
% --------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{DissertateUSU}[2017/07/01 v1.0 DissertateUSU Class]
\LoadClass[12pt, oneside, letterpaper]{book}

%%%%%%%%%%%%%%%%%
%%   Options   %%
%%%%%%%%%%%%%%%%%
\RequirePackage{etoolbox}

%{{{ Line-spacing control
\RequirePackage{setspace}
\newcommand{\singlespaced}{\singlespacing}
\newcommand{\doublespaced}{\doublespacing}
\newcommand{\triplespaced}{\setstretch {2.6}}

%}}}

% Text layout.
\usepackage{ragged2e}
\RaggedRight
\RequirePackage{graphicx}
\usepackage{fixltx2e}
\usepackage{parskip}
\setlength{\parindent}{1cm}
\RequirePackage{lettrine}
\RequirePackage{setspace}
\RequirePackage{verbatim}
\usepackage{indentfirst}
\usepackage{bibentry}


% Fonts.
\RequirePackage{color}
\RequirePackage{xcolor}
\usepackage{hyperref}
\RequirePackage{url}
\RequirePackage{amssymb}
\RequirePackage{mathspec}

\AtBeginEnvironment{tabular}{\addfontfeature{RawFeature=+tnum}}
\widowpenalty=300
\clubpenalty=300
% \setromanfont[Numbers=OldStyle, Ligatures={Common, TeX}, Scale=1.0]{EB Garamond}
% \newfontfamily{\smallcaps}[RawFeature={+c2sc,+scmp}]{EB Garamond}
% \setsansfont[Scale=MatchLowercase, BoldFont={Lato Bold}]{Lato Regular}
% \setmonofont[Scale=MatchLowercase]{Source Code Pro}
\RequirePackage[labelfont={bf,small,singlespacing},
                textfont={small,singlespacing},
                justification={justified,RaggedRight},
                singlelinecheck=false,
                margin=0pt,
                figurewithin=chapter,
                tablewithin=chapter]{caption}
\RequirePackage{microtype}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Headings & Headers         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{fancyhdr}
\RequirePackage[tiny, md, sc]{titlesec}
\setlength{\headheight}{15pt}
\RequirePackage{titling}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Page number in lower center.     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\fancypagestyle{plain}{%
\fancyhf{} % clear all header and footer fields
\fancyfoot[C]{\thepage} % except the right
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAYOUT: TABLE of CONTENTS -auto fill %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setcounter{tocdepth}{2}
\usepackage[titles]{tocloft}
\usepackage[titletoc]{appendix}
\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}
\renewcommand{\cftchapfont}{\normalsize \normalfont}
\renewcommand\contentsname{CONTENTS}
\renewcommand\listfigurename{LIST OF FIGURES}
\renewcommand\listtablename{LIST OF TABLES}
\renewcommand{\cftchapafterpnum}{\vspace{10pt}}
\newcommand{\hsp}{\hspace{20pt}}

\def\boldcmd{\bfseries}

\def\typesetHeading#1{{\normalfont\large\boldcmd{#1}}}
\def\typesetColumnHeading#1{{\boldcmd #1}}
\def\typesetChapterTitle#1{{#1}}

\newcommand\vitaname{Vita}
\newcommand\vitanameLC{Vita}
\newcommand\publicationsname{Publications}
\newcommand\fieldsname{Fields of Study}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%.        VITA                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{vita}{
  \chapter*{\typesetHeading\vitaname
        \@mkboth{\uppercase{\vitaname}}{\uppercase{\vitaname}}}%
  \addcontentsline{toc}{section}{{\vitanameLC}}
  \thispagestyle{plain}
  \pagestyle{plain}
  \vskip 10pt
  \singlespaced
}
{\par\doublespaced}
%%
%%  \dateitem{date}{what}
%%
%%  Creates a listing of important dates. The date and what are 
%%  separated by a \dotfill and the what is in a parbox
%%
\def\dateitem#1#2{\noindent #1 \dotfill \parbox[t]{2.75in}{#2}\par~\par}
%%
%%  \begin{publist} ...\end{publist}
%%
%%  Used around the list of publications.  Items done and separated
%%  using the commands below
%%
\newenvironment{publist}%
  {\vspace*{15pt}\begin{center}{\typesetHeading{\publicationsname}}%
   \end{center}\vspace*{10pt}}
  {\vspace{3pt}}
%%
%% \researchpubs  and \instructpubs
%%
%% Used to create ``sub-headings'' within the list of publications
%%
\def\researchpubs{\noindent {\bf Research Publications}\par~~\par}
\def\instructpubs{\noindent {\bf Instructional Publications}\par~~\par}
%%
%%  \pubitem{reference}
%%
%%  Used to list the references you have within a publist environment
%%  Note:  Easiest to use bibtex to generate this list.
%%
\def\pubitem#1{\noindent #1 \par~~\par}

%%
%% \begin{fieldsstudy} ...\end{fieldsstudy}
%%
%%  Used to create the list of your major and what fields you have 
%%  looked at within your major.
%%
\newenvironment{fieldsstudy}%
  {\begin{center}{\typesetHeading\fieldsname}\end{center}%
   \vspace*{15pt}}
  {}
% ted - (2010) default to \@unit if star present
%%
%% \majorfield{Field}
%% \majorfield*
%%
%%  Used to create a line indicating your major
%%    e.g. Computer Science and Engineering
%%
\def\unit#1{\def\@unit{#1}}
\def\@unit{Pharmaceutical Sciences}
\def\graduationyear#1{\def\@graduationyear{#1}}
\def\@graduationyear{\year}
\def\majorfield{\@ifstar{\@majorfield{\@unit}}{\@majorfield}}
\def\@majorfield#1{\noindent Major Field:  #1\par~\par}
%%
%% \onestudy{Topic}{Professor}
%%
%%  Used to create a single line if you have studies only one topic.
%%  If have looked at more than one thing, should use the environment below.
%%
\def\onestudy#1#2{Studies in #1: #2}
%%
%%  The following is used to list the topics you have studies and who you
%%  studied with.
%%
%% \begin{studieslist}
%% \studyitem{Topic}{Professor}
%% \studyitem{Topic}{Professor}
%%    ...
%% \studyitem{Topic}{Professor}
%% \end{studieslist}
%%

\def\studieslist{Studies in:\par\vspace*{8pt}\hspace*{2em}\begin{tabular}{l l}}
\def\endstudieslist{\end{tabular}}

\def\studyitem#1#2{#1 & #2\\}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAYOUT: List of TABLES - auto fill %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{{{ List of Tables
% LoT should be single-spaced within entries and double-spaced between entries
\renewcommand*{\listoftables}{
    \singlespaced
    % \vskip0pt
    \chapter*{ \begin{center} \textnormal{LIST OF TABLES}  \end{center}}
    \addcontentsline{toc}{section}{LIST OF TABLES}
    % \addcontentsline{toc}{chapter}{\textnormal{LIST OF TABLES}}
    % \addvspace{10pt}
    % \typeout{!! List of Tables...!!}
    \noindent Table \hfill Page\\         % newline allows us to use a vspace...
    % vspace adjusts spacing to get a doublespace between table headings
    % and table contents
    % \vspace{\baselineskip}
    \@starttoc{lot}     
    \newpage
    % \doublespaced
}
% Ensure double-spacing between entries
\renewcommand*\l@table{\par\addvspace{\bigskipamount}\pagebreak[2]
 \@dottedtocline{1}{1.5em}{2.3em}}
%}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAYOUT: List of FIGURES -auto fill %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{{{ List of Figures
% LoF should be single-spaced within entries and double-spaced between entries 
\renewcommand*{\listoffigures}{
    \singlespaced
    % \vskip0pt
    \chapter*{\begin{center} \textnormal{LIST OF FIGURES} \end{center}}
    \addcontentsline{toc}{section}{LIST OF FIGURES}
    % \addcontentsline{toc}{chapter}{\textnormal{LIST OF FIGURES}}
    % \addvspace{10pt} 
    % \typeout{!! List of Figures...!!}
    \noindent Figure \hfill Page\\         % newline allows us to use a vspace...
    % vspace adjusts spacing to get a doublespace between table headings
    % and table contents
    % \vspace{\baselineskip} 
    \@starttoc{lof}     
    \newpage 
    % \doublespaced
}
% Ensure double-spacing between entries
\renewcommand*\l@figure{\par\addvspace{\bigskipamount}\pagebreak[2]
 \@dottedtocline{1}{1.5em}{2.3em}}
%}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAYOUT: SECTION HEADINGS- auto fill %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\setthesection}{\roman{chapter}.A\roman{section}}
\usepackage{quotchap}
\parskip 0pt

\usepackage{titlesec, blindtext, color}

\titleformat{\chapter}[hang]
{\vspace{-60pt}\center}   %% spacing of chapter headings
{{CHAPTER \thechapter:}}
{1em}
{}


\titleformat{\section}[display]
{\normalfont\normalsize\bfseries}
{}
{0pt}
{{}}

\titleformat{\subsection}[display]
{\normalfont\normalsize\itshape}
{}
{0pt}
{{}}

\titleformat{\subsubsection}[display]
{\normalfont\normalsize}
{}
{0pt}
{{}}

% An environment for paragraph-style section.
\providecommand\newthought[1]{%
   \addvspace{1.0\baselineskip plus 0.5ex minus 0.2ex}%
   \textsc{#1}%
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure numbering - chapter.number %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\thefigure}{\arabic{section}.\arabic{figure}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAYOUT: Title Page - info filled above    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\maketitle}{
    \thispagestyle{empty}
    \vspace*{35pt}
    \begin{center}
    \doublespaced
    \MakeUppercase{\@title}\\
    \vspace{35pt}
    Dissertation \\
    \vspace{35pt} 
    Presented in Partial Fulfillment of the Requirements for the Degree Doctor of Philosophy in the Graduate School of The Ohio State University \\
    \vspace{20pt}
    By \\
    \singlespaced
    \@author \\
    \singlespaced
    Graduate Program in\degree \\
    \vspace{35pt} 
    \gradschoolguy \\
    \singlespaced
    \yeardegree  \\ 
    \vspace{46pt}
    Dissertation Committee: \\
    \singlespaced
    \committeeone \\
    \singlespaced
    \committeetwo \\
    \singlespaced
    \committeethree \\
    \singlespaced
    \committeefour \\
    \end{center}

    \clearpage
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAYOUT: Copy Right - info filled above    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\copyrightpage}{
	\vspace*{\fill}
  \begin{center}
	\doublespacing
	Copyright \hspace{3pt} 
	  \scshape \small \copyright  \hspace{3pt} 
	  Abbas A Rizvi \hspace{3pt}\yeardegree\\
	All Rights Reserved
  \end{center}
	\vspace*{\fill}
}


