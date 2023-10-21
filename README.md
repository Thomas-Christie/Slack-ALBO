# Slack-ALBO

Repository for running experiments from [Picheny et al.
2016](https://arxiv.org/pdf/1605.09466.pdf). The code used to run the
experiments presented in our paper (on the `ackley_10` and `lockwood` problems), and the
results, can be found in the `final_results/ackley_10/paper_results` and
`final_results/lockwood/paper_results` directories. Within these directories, the
results are in the `data` directories and the source code used to obtain them is
provided in the `details.md` file. Source code is also found in the `ackley10.R` and
`lockwood/runlock/lockwood_al_runner.R` files. Note that we made some edits to the
original code provided in the [laGP package](https://bobby.gramacy.com/r_packages/laGP/)
in order to fix several (minor) bugs we found in the code. 

In order to run the `Lockwood` problem, you will need to build the executable for
running the simulator using the code provided on Robert Gramacy's website
[here](https://bobby.gramacy.com/surrogates/). You will then need to utilise the
`lockwood_al_runner.R` file provided in the `lockwood/runlock/` directory of this repository. 
