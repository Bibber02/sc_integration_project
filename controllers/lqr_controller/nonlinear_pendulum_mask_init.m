classdef nonlinear_pendulum_mask_init

    methods(Static)

        % Following properties of 'maskInitContext' are available to use:
        %  - BlockHandle
        %  - MaskObject
        %  - MaskWorkspace: Use get/set APIs to work with mask workspace.
        function MaskInitialization(maskInitContext)
            % ---------------------------------------------------------
            % Mask initialization for the rotational pendulum Plant
            % subsystem. Runs every time the mask is initialized.
            %
            % Steps:
            %   1. Verify that the identified parameters are present in
            %      the base workspace (set by setup_pendulum_model.m).
            %   2. Read the four initial-state mask parameters.
            %   3. Promote them into the mask workspace under the names
            %      that the Integrator blocks reference.
            % ---------------------------------------------------------

            % --- 1. Check base workspace for required parameters ----
            requiredVars = {'p_alpha','p_b1','p_c1','p_g1','p_u','p_0', ...
                            'p_b2','p_g2','p_c2','p_s_d2','v_s2', ...
                            'eps_v1','eps_v2','l_1','g','p_c'};

            missing = {};
            for k = 1:numel(requiredVars)
                if ~evalin('base', sprintf('exist(''%s'',''var'')', requiredVars{k}))
                    missing{end+1} = requiredVars{k};
                end
            end

            if ~isempty(missing)
                error(['Missing workspace variables: %s\n' ...
                       'Run setup_pendulum_model.m before simulating.'], ...
                       strjoin(missing, ', '));
            end

            % --- 2. Read mask parameters -----------------------------
            mw = maskInitContext.MaskWorkspace;
            theta_1_0     = mw.get('theta_1_0_mask');
            theta_2_0     = mw.get('theta_2_0_mask');
            theta_1_dot_0 = mw.get('theta_1_dot_0_mask');
            theta_2_dot_0 = mw.get('theta_2_dot_0_mask');

            % --- 3. Write the resolved IC names into the mask ws ----
            % Integrator blocks reference these names directly.
            mw.set('theta_1_0',     theta_1_0);
            mw.set('theta_2_0',     theta_2_0);
            mw.set('theta_1_dot_0', theta_1_dot_0);
            mw.set('theta_2_dot_0', theta_2_dot_0);
        end

        % Per-parameter callbacks fire when the user edits that
        % specific dialog field. Empty by default; populate to add
        % per-field validation if desired.

        function theta_1_0_mask(callbackContext) %#ok<INUSD>

        end

        function theta_2_0_mask(callbackContext) %#ok<INUSD>

        end

        function theta_1_dot_0_mask(callbackContext) %#ok<INUSD>

        end

        function theta_2_dot_0_mask(callbackContext) %#ok<INUSD>

        end
    end
end