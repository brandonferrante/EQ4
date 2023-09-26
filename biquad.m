classdef biquad < audioPlugin
    properties
        HS_FREQ = 10000; %high shelf frequency
        HS_GAIN = 0;     %high shelf gain

        HMF_FREQ = 5000;    %hi-mid frequency
        HMF_GAIN = 0;       %hi-mid gain

        LMF_FREQ = 500;     %low-mid frequency
        LMF_GAIN = 0;       %low-mid gain

        HPF_FREQ = 50;      %high pass frequency

        fs = 44100;         %sample rate
        fn = 22050;          %Nyquist

        BYPASS = 'off';     %bypass ON or OFF

    end

    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('HS_FREQ', 'DisplayName', ...
                'Hi-Shelf Frequency','Label', 'Hz', 'Mapping', ...
                 {'log', 2500, 20000}),...
            audioPluginParameter('HS_GAIN', 'DisplayName', ...
                'Hi-Shelf Gain', 'Label', 'dB', 'Mapping', ...
                 {'lin', -12, 12}),...
            audioPluginParameter('HMF_FREQ', 'DisplayName', ...
                'Hi-Mid Frequency', 'Label', 'Hz', 'Mapping', ...
                 {'log', 800, 12500}),...
            audioPluginParameter('HMF_GAIN', 'DisplayName', ...
                'Hi-Mid Gain', 'Label', 'dB', 'Mapping', ...
                 {'lin', -12, 12}),...
            audioPluginParameter('LMF_FREQ', 'DisplayName', ...
                'Lo-Mid Frequency', 'Label', 'Hz', 'Mapping', ...
                 {'log', 75, 1000}),...                 
            audioPluginParameter('LMF_GAIN', 'DisplayName', ...
                'Lo-Mid Gain', 'Label', 'dB', 'Mapping', ...
                 {'lin', -12, 12}),...
            audioPluginParameter('HPF_FREQ', 'DisplayName', ...
                'High Pass', 'Label', 'Hz', 'Mapping', ...
                 {'log', 30, 400}),...      
            audioPluginParameter('BYPASS', 'DisplayName', ...
                'Bypass','Mapping', ...
                 {'enum', 'off', 'on'}));                
    end

    properties (Access = private)
        filter_HS = struct('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_HMF = struct('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_LMF = struct('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_HPF = struct('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);

    end

    methods

        function out = process(plugin, in)
            
            out = zeros(size(in));

            for ch = 1:min(size(in))

                x = in(:,ch);

                [y1, plugin.filter_HS.w(:,ch)] = processBiquad(x, plugin.filter_HS, ch);
                [y2, plugin.filter_HMF.w(:,ch)] = processBiquad(y1, plugin.filter_HMF, ch);
                [y3, plugin.filter_LMF.w(:,ch)] = processBiquad(y2, plugin.filter_LMF, ch);
                [y4, plugin.filter_HPF.w(:,ch)] = processBiquad(y3, plugin.filter_HPF, ch);

                if strcmp(plugin.BYPASS, 'on')
                    out(:,ch) = x;

                else
                    out(:,ch) = y4;

                end

            end

        end

        % DAW reset
        function reset(plugin)

            plugin.fs = getSampleRate(plugin);
            plugin.fn = plugin.fs/2;

            plugin.filter_HS.w = [0 0; 0 0];
            plugin.filter_HMF.w = [0 0; 0 0];
            plugin.filter_LMF.w = [0 0; 0 0];
            plugin.filter_HPF.w = [0 0; 0 0];

        end


%Hi shelf

        % set Hi-shelf freq
        function set.HS_FREQ(plugin, val)

            plugin.HS_FREQ = val;
            update_HS(plugin);

        end

        % set Hi-shelf gain
        function set.HS_GAIN(plugin, val)
            
            plugin.HS_GAIN = val;
            update_HS(plugin);

        end

        function update_HS(plugin)

            Q = 0.5;
            f0 = plugin.HS_FREQ;
            gain = plugin.HS_GAIN;
            w0 = 2 * pi * f0 / plugin.fs; %convert to radians/sec
            alpha = sin(w0)/2*Q;
            A = sqrt(db2mag(gain));

            %HS coefficient calculator 
            plugin.filter_HS.a0 =    A*( (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha );
            plugin.filter_HS.a1 = -2*A*( (A-1) + (A+1)*cos(w0));
            plugin.filter_HS.a2 =    A*( (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha );
            plugin.filter_HS.b0 =        (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha;
            plugin.filter_HS.b1 =    2*( (A-1) - (A+1)*cos(w0));
            plugin.filter_HS.b2 =        (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha;

        end

%Hi-mid

        % set Hi-mid freq
        function set.HMF_FREQ(plugin, val)

            plugin.HMF_FREQ = val;
            update_HMF(plugin);

        end

        % set Hi-mid gain
        function set.HMF_GAIN(plugin, val)
            
            plugin.HMF_GAIN = val;
            update_HMF(plugin);

        end

        function update_HMF(plugin)

            Q = 0.5;
            f0 = plugin.HMF_FREQ;
            gain = plugin.HMF_GAIN;
            w0 = 2 * pi * f0 / plugin.fs; %convert to radians/sec
            alpha = sin(w0)/2*Q;
            A = sqrt(db2mag(gain));

            % HMF coefficient calculator 
            plugin.filter_HMF.a0 =  1 + alpha*A;
            plugin.filter_HMF.a1 = -2*cos(w0);
            plugin.filter_HMF.a2 =  1 - alpha*A;
            plugin.filter_HMF.b0 =  1 + alpha/A;
            plugin.filter_HMF.b1 = -2*cos(w0);
            plugin.filter_HMF.b2 =  1 - alpha/A;

        end

%Low-mid

        % set Low-mid freq
        function set.LMF_FREQ(plugin, val)

            plugin.LMF_FREQ = val;
            update_LMF(plugin);

        end

        % set Lo-mid gain
        function set.LMF_GAIN(plugin, val)
            
            plugin.LMF_GAIN = val;
            update_LMF(plugin);

        end

        function update_LMF(plugin)

            Q = 0.5;
            f0 = plugin.LMF_FREQ;
            gain = plugin.LMF_GAIN;
            w0 = 2 * pi * f0 / plugin.fs; %convert to radians/sec
            alpha = sin(w0)/2*Q;
            A = sqrt(db2mag(gain));

            %LMF coefficient calculator 
            plugin.filter_LMF.a0 =  1 + alpha*A;
            plugin.filter_LMF.a1 = -2*cos(w0);
            plugin.filter_LMF.a2 =  1 - alpha*A;
            plugin.filter_LMF.b0 =  1 + alpha/A;
            plugin.filter_LMF.b1 = -2*cos(w0);
            plugin.filter_LMF.b2 =  1 - alpha/A;

        end
        
%High pass

        % set Hi-pass freq
        function set.HPF_FREQ(plugin, val)

            plugin.HPF_FREQ = val;
            update_HPF(plugin);

        end

        function update_HPF(plugin)

            Q = 0.5;
            f0 = plugin.HPF_FREQ;
            w0 = 2 * pi * f0 / plugin.fs; %convert to radians/sec
            alpha = sin(w0)/2*Q;

            %HPF coefficient calculator 
            plugin.filter_HPF.a0 =  (1 + cos(w0))/2;
            plugin.filter_HPF.a1 = -(1 + cos(w0));
            plugin.filter_HPF.a2 =  (1 + cos(w0))/2;
            plugin.filter_HPF.b0 =   1 + alpha;
            plugin.filter_HPF.b1 =  -2*cos(w0);
            plugin.filter_HPF.b2 =   1 - alpha;

        end

        %set bypass
        function set.BYPASS(plugin, val)

            plugin.BYPASS = val;

        end

    end

end