% Copyright (c) 2026 Ziyi Yang (ziyi016@e.ntu.edu.sg).
% Extended from the MAML/FxLMS code lineage of Shi Dongyuan's Meta repository:
% https://github.com/ShiDongyuan/Meta
%
classdef MAML_Nstep_forget_S2
    properties
        Psi   % The initial secondary-path filter (meta init)
    end
    methods
        function obj = MAML_Nstep_forget_S2(len_s)
            % len_s : the length of secondary-path FIR
            obj.Psi = zeros(len_s,1);
        end

        function [obj, S_task, ErS] = MAML_initial_S(obj, U, Y, muS, lamdaS, epslonS)
            % U : identification input (excitation) segment (length = len_s)
            % Y : identification output (measured at error mic) segment
            % muS, lamdaS, epslonS: step size, forgetting factor, meta rate.
            U   = flipud(U);
            Y   = flipud(Y);
            GradS = 0;             % temporal gradient accumulator
            ErS   = 0;
            Ls    = length(obj.Psi);

            % <--A1--> one inner-step from meta init
            e  = Y(1) - obj.Psi' * U;
            S_task = obj.Psi + muS * e * U;  % assumed task-optimal after one step

            % <--A2--> accumulate meta gradient with forgetting
            for jj = 1:Ls
                Ud = [U(jj:end); zeros(jj-1,1)];
                e  = Y(jj) - S_task' * Ud;
                GradS = GradS + epslonS * (muS/Ls) * e * Ud * (lamdaS^(jj-1));
                if jj == 1, ErS = e; end
            end

            % <--A3--> meta update of Psi
            obj.Psi = obj.Psi + GradS;
        end
    end
end
