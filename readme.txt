optimal power flow microgrids test feeder
Edits from the SDP model: (1) because of the transformer ratio multiple solution problem, the transformer ratio in non-convex model solved by GAMS is equal to the results calculated by SDP model through MATLAB.
The code requires (1) MATLAB,GAMS (2) CVX to call the solver (3) SDPT3 as the SDP solver.
Reference: Three-Phase Optimal Power Flow for Networked Microgrids Based on Semi-definite Programming Convex Relaxation