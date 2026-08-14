#!/usr/bin/env bash

# Exit immediately if a command fails
set -e

PROJECT_DIR="BLM1995"

echo "🚀 Creating modular LaTeX directory structure for '$PROJECT_DIR'..."

# Create directory tree
mkdir -p "$PROJECT_DIR/sections"
mkdir -p "$PROJECT_DIR/tables"
mkdir -p "$PROJECT_DIR/figures"

# ---------------------------------------------------------------------
# 1. Headers File (blm_headers.tex)
# ---------------------------------------------------------------------
cat << 'EOF' > "$PROJECT_DIR/blm_headers.tex"
% =====================================================================
% Document Class & Required Packages
% =====================================================================
\documentclass[11pt,a4paper]{article}

\usepackage[utf8]{utf8}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{bm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{cite}
\usepackage{microtype}
\usepackage{hyperref}

\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    citecolor=blue,
    urlcolor=blue
}
EOF

# ---------------------------------------------------------------------
# 2. Macros File (blm_macros.tex)
# ---------------------------------------------------------------------
cat << 'EOF' > "$PROJECT_DIR/blm_macros.tex"
% =====================================================================
% Boundary-Layer Meteorology Macros
% =====================================================================

% Stability Parameters
\newcommand{\Ri}{\mathrm{Ri}}
\newcommand{\Ric}{\mathrm{Ri}_c}
\newcommand{\Rib}{\mathrm{Ri}_b}
\newcommand{\Rig}{\mathrm{Ri}_g}
\newcommand{\Rf}{\mathrm{R}_f}
\newcommand{\atheta}{\alpha_\theta}
\newcommand{\Ria}{\left(\frac{\Ri}{\atheta}\right)}
\newcommand{\Riasq}{\left(\frac{\Ri}{\atheta}\right)^2}

% Stability Functions
\newcommand{\fm}{f_m}
\newcommand{\fh}{f_h}
\newcommand{\sfm}{\sqrt{f_m}}
\newcommand{\sfh}{\sqrt{f_h}}
\newcommand{\Psim}{\Psi_m}
\newcommand{\Psih}{\Psi_h}

% Dimensionless Gradient Functions
\newcommand{\phim}{\phi_m}
\newcommand{\phih}{\phi_h}

% Layer Coefficients
\newcommand{\betam}{\beta_m}
\newcommand{\betah}{\beta_h}
\newcommand{\ubm}{b_m}
\newcommand{\ubh}{b_h}

% Prandtl Numbers & Exchange Coefficients
\newcommand{\Prt}{\mathrm{Pr}_t}
\newcommand{\Przero}{\mathrm{Pr}_0}
\newcommand{\Km}{K_m}
\newcommand{\Kh}{K_h}

% Turbulence & Length Scales
\newcommand{\ustar}{u_*}
\newcommand{\thetastar}{\theta_*}
\newcommand{\wstar}{w_*}
\newcommand{\vk}{\kappa}
\newcommand{\zo}{z_0}
\newcommand{\zone}{z_1}
\newcommand{\zbar}{\bar{z}}
\newcommand{\LMO}{L}

% Monin--Obukhov & Regimes
\newcommand{\MO}{Monin--Obukhov}
\newcommand{\MOST}{Monin--Obukhov Similarity Theory}
\newcommand{\zetaMO}{\zeta}
\newcommand{\SBL}{Stable Boundary Layer}
\newcommand{\UBL}{Unstable Boundary Layer}

% Temperature & Vectors
\newcommand{\thetav}{\theta_v}
\newcommand{\thetag}{\theta_g}
\newcommand{\Vone}{\mathbf{V}_1}
\newcommand{\Vonesq}{|\mathbf{V}_1|^2}

% Repeated Fractions
\newcommand{\lnz}{\ln\left(\frac{\zone}{\zo}\right)}
\newcommand{\lnzsq}{\left(\ln\frac{\zone}{\zo}\right)^2}

% Math Operators
\newcommand{\dd}{\,\mathrm{d}}
\newcommand{\pd}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\od}[2]{\frac{\mathrm{d}#1}{\mathrm{d}#2}}
\newcommand{\vect}[1]{\bm{#1}}
\newcommand{\grad}{\nabla}
\newcommand{\Div}{\nabla\!\cdot}
\newcommand{\Curl}{\nabla\times}
EOF

# ---------------------------------------------------------------------
# 3. Main Document (main.tex)
# ---------------------------------------------------------------------
cat << 'EOF' > "$PROJECT_DIR/main.tex"
\input{blm_headers}
\input{blm_macros}

\title{Analytical Solutions for Boundary-Layer Stability Functions}
\author{Boundary-Layer Research Group}
\date{\today}

\begin{document}

\maketitle

\begin{abstract}
Abstract text goes here...
\end{abstract}

\input{sections/introduction}
\input{sections/exchange_coefficients}
\input{sections/stability_functions}
\input{sections/discussion}
\input{sections/conclusions}

\appendix
\input{sections/appendix}

\bibliographystyle{plain}
\bibliography{references}

\end{document}
EOF

# ---------------------------------------------------------------------
# 4. Section Placeholders
# ---------------------------------------------------------------------
cat << 'EOF' > "$PROJECT_DIR/sections/introduction.tex"
\section{Introduction}
Introduction section content goes here...
EOF

cat << 'EOF' > "$PROJECT_DIR/sections/exchange_coefficients.tex"
\section{Exchange Coefficients}
Exchange coefficients content goes here...
EOF

cat << 'EOF' > "$PROJECT_DIR/sections/stability_functions.tex"
\section{Stability Functions}
Stability functions content goes here...
EOF

cat << 'EOF' > "$PROJECT_DIR/sections/discussion.tex"
\section{Discussion}
Discussion content goes here...
EOF

cat << 'EOF' > "$PROJECT_DIR/sections/conclusions.tex"
\section{Conclusions}
Conclusions content goes here...
EOF

cat << 'EOF' > "$PROJECT_DIR/sections/appendix.tex"
\section*{Appendix: Analytical Solution of Unstable Cubic}
Appendix content goes here...
EOF

# ---------------------------------------------------------------------
# 5. Tables & References Placeholders
# ---------------------------------------------------------------------
cat << 'EOF' > "$PROJECT_DIR/tables/table1.tex"
% Table 1 Placeholder
\begin{table}[ht]
\centering
\caption{Boundary layer parameter comparison.}
\begin{tabular}{lll}
\toprule
Parameter & Description & Default Value \\
\midrule
$\ubm$ & Unstable momentum coefficient & 15.0 \\
$\ubh$ & Unstable heat coefficient & 9.0 \\
\bottomrule
\end{tabular}
\end{table}
EOF

cat << 'EOF' > "$PROJECT_DIR/references.bib"
@article{MoninObukhov1954,
  author  = {Monin, A. S. and Obukhov, A. M.},
  title   = {Basic laws of turbulent mixing in the atmosphere near the ground},
  journal = {Tr. Akad. Nauk SSSR Geofiz. Inst.},
  volume  = {24},
  pages   = {163--187},
  year    = {1954}
}
EOF

echo "✨ Done! Project set up successfully under './$PROJECT_DIR'."