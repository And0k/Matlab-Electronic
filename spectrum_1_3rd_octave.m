%% Spectrum of sound pressure level (SPL) from Oscilloscope Data Files
%% Copyright © 2022 Andrey Korzh <ao.korzh@gmail.com>
global save_dir lim_x lim_y

%% User configuration
dir_data = 'D:\WorkData\_experiment\underwater_noise\_raw';

% Spectrum calculation and drawing limits
lim_x = [3 10000];  % [min, max] Hz, if max > 10000 then still 10000 will be used
lim_y = [50 110];  % dB ref: 20 uPa

% Amplification coefficients applied to a sound signal
coef_hydrophone = 0.53;           % pC/Pa
coef_preamplifier = 100;          % mV/pC (use 1, 10 or 100)

%% Start
% Coefficient to convert recorded voltage to a sound pressure in Pa
coef_V_to_Pa = 1 / (coef_hydrophone * coef_preamplifier * 1E-3);  % Pa/V

b_save = true;  % if false then no files will be saved

[dir_parent, ~, ~] = fileparts(dir_data);
save_dir = [dir_parent '/data_out'];
if ~exist(save_dir, "dir")
    mkdir(save_dir);
end

% Unzip all zip files in same folder (dir_data), delete zip
file_names = dir([dir_data '/*.zip']);
if ~isempty(file_names)
    fprintf(1, "Extracting files from %d zip archive(s).\n", length(file_names))
    for file_name_struct = file_names'
        file_name = file_name_struct.name;
        unzip([dir_data '/' file_name], dir_data);
        delete([dir_data '/' file_name]);
    end
end

% Process all found *.txt files
file_names = dir([dir_data '/*.txt']);
fprintf(1, "Processing %d found file(s)", length(file_names))
i = 0;
for file_name_struct = file_names'
    file_name = file_name_struct.name;
    i = i + 1;
    fprintf(1, '\n%d. ', i)
    [ar, fs, t_start] = read_data([dir_data '/' file_name], coef_V_to_Pa);
    data_title = sprintf('%s (%d samples recorded %s)', file_name, length(ar), t_start);
    fprintf(1, data_title)

    if false  % set true to display raw data
        figure('Name', ['Raw data ' data_title]) %#ok<UNRCH> 
        plot(ar, '-', 'DisplayName', 'raw');
        % legend(ax, 'Box','off', 'Color','none')
        % legend('show')
    end

    get_spectrum(ar, fs, file_name, data_title, b_save);
end
fprintf(1, "Ok");

%% Functions

function get_spectrum(ar, fs, filename, data_title, b_save)
    % Calculate spectrum
    % fs: sampling frequency
    % data_title: figure title end. If no then do not create figure, required if want save image
    % b_save: whether to save result (figure image and csv)
    global save_dir lim_x lim_y     % lim_x maximum is set here to be <= 10000 Hz
    persistent fs_prev filter_bank  % filter bank to calculate high resolution spectrum depends on fs
    persistent ha_octave ha_psd     % handles of axes to reuse
    pref = 0.00002;     % 20 micro Pascals reference sound pressure level (0 dB SPL = 0.00002 Pa)
    SPL_units_label = ', dB re 20 uPa';
    SPL_units_label_utf8 = ', dB re 20 μPa';

    % units of PSD (not for octave spectrum which we hard-coded to be Hz in poctave_localplot_10lg())
    x_units_plot = 'kHz';
    x_units_plot_coef = 1E-3;

    if isnan(lim_x(end)) || lim_x(end) > 10000
        lim_x(end) = 10000;  % Hz
    end

    % Spectrogram
    %%%%%%%%%%%%%
    if false  % data_title
        figure('Name', ['Spectrogram spectrum of ' data_title]) %#ok<UNRCH> 
        poctave(ar / pref, fs, 'spectrogram', 'OverlapPercent', 50)
        title([get(gca, 'title').String '. Input: ' filename])
        if b_save
            path_spectrogram = [save_dir '/spectrogram_' filename(1 : end-4) '.jpg'];
            print(gcf, '-djpeg', path_spectrogram, '-noui');
        end
    end

    % Octave spectrum
    %%%%%%%%%%%%%%%%%
    BandsPerOctave = 3;
    % Average power over (sub)octave bands and center frequencies
    [Power, freq_oct] = poctave(ar, fs, ...
        'BandsPerOctave', BandsPerOctave, ...
        'FrequencyLimits', lim_x);
    % PSD
    freq_edges = [freq_oct.*2^(-1/(2*BandsPerOctave)); freq_oct(end)*2^(1/(2*BandsPerOctave))];
    bands_oct = diff(freq_edges);
    psd_oct = Power ./ bands_oct;

    if data_title
        if isempty(ha_octave) || ~ishandle(ha_octave)
            hf = figure('Name', ['1/3-octave power spectrum of ' data_title]);      
            grid on  % also creates axes on figure
            ha_octave = get(hf, 'CurrentAxes');
            set(ha_octave, 'NextPlot','replacechildren');
            ylabel(ha_octave, ['Sound Pressure Level' SPL_units_label_utf8]);
            % x label will be added by poctave_localplot_10lg()
        else
            hf = ha_octave.Parent;
            set(hf, 'Name', ['1/3-octave power spectrum of ' data_title]);
        end
        poctave_localplot_10lg(ha_octave, Power / pref.^2, freq_oct);
        % poctave(ar / pref, fs, 'BandsPerOctave', 3, 'FrequencyLimits',[max(3, lim_x(1)), lim_x(2)]);
        title(['1/3-octave spectrum. ' filename], 'Parent', ha_octave);
        if b_save
            path_name_psd = [save_dir '/spectr_oct_one3rd_' filename(1 : end-4)];
            print(hf, '-djpeg', [path_name_psd, '.jpg'], '-noui');
        end
    end
    
    % Save to text file
    if b_save
        % Convert to calibrated SPL values (referenced to 20 uPa)
        % see also https://www.mathworks.com/matlabcentral/answers/431461-poctave-return-value-for-acoustics-analysis
        SPL = 10 * log10(Power / pref.^2);
        writetable( ...
            table( ...
                freq_oct, SPL, 'VariableNames', { ...  % bands_oct, psd_oct, 
                    'freq, Hz',  ...  % 'band, Hz', 'PSD, Pa^2/Hz',
                    ['SPL' SPL_units_label] ...
                    } ...
                ), ...
            [path_name_psd '.csv'] ...
            )
    end

    % High resolution spectrum with filter bank
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(fs_prev) || fs ~= fs_prev
        fs_prev = fs;
        filter_bank = get_spectrum_analyzer(fs);
        if isempty(filter_bank)  % if failed
            fs_prev = NaN;       % can not use previous filter_bank later
            return
        end
    end
    filter_bank(ar / pref);
    % set(gca().Color, [1,1,1]) How???
    while ~filter_bank.isNewDataReady
        drawnow
    end
    release(filter_bank)  % displays output in live editor temporary?
    tbl_data = getSpectrumData(filter_bank);

    % PSD
    bands = diff(tbl_data.FrequencyVector{1}(1:2));  % all the same
    psd = tbl_data.Spectrum{1} .^ 2 / bands; % Power = tbl_data.Spectrum{1} .^ 2;

    % noise dBm:  10*log10((noiseVar/(NumFreqBins/2))/1e-3)
    % sqrt((10.^(tbl_data.Spectrum{1}/10))*1e-3); Vrms from dBm:
    % SPL = tbl_data.Spectrum{1} + 10 * log10(1E-3/(pref^2));  % dBm to dB SPL
    SPL = 20 * log10(tbl_data.Spectrum{1} / pref);  % Vrms to SPL
    % SPL = 10 * log10(tbl_data.Spectrum{1} / pref^2);  % Power to SPL
    if data_title
        if isempty(ha_psd) || ~ishandle(ha_psd)
            hf = figure('Name', ['High resolution spectrum of ' data_title]);
            ha_psd = newplot();
            set(ha_psd, 'NextPlot','replacechildren', ...
                'YScale', 'log');
            xlabel(ha_psd, ['Frequency, ' x_units_plot]);
            ylabel(ha_psd, 'Power Spectrum Density, Pa^2/Hz');
            xlim(lim_x * x_units_plot_coef);  % kHz
            % ylim(lim_y);
            grid(ha_psd, 'on');
            legend(ha_psd, 'Box','off', 'Color','none')
            legend(ha_psd, 'show')
        else
            hf = ha_psd.Parent;
            set(hf, 'Name', ['High resolution spectrum of ' data_title]);
        end
        semilogx(ha_psd, ...
            tbl_data.FrequencyVector{1} * x_units_plot_coef, psd*1E-9, '-', ...
            'DisplayName', 'uniform intervals', ...
            'Color', 'k' ...
            );
        hold(ha_psd, 'on');
        plot(ha_psd, freq_oct * x_units_plot_coef, psd_oct, '-', ...
            'DisplayName', '1/3 octave intervals', ...
            'Color', 'r' ...
            );
        set(ha_psd, 'NextPlot','replacechildren');
        title(['High resolution spectrum. ' filename], 'Parent', ha_psd);
        if b_save
            path_name_psd = [save_dir '/spectr_hr_' filename(1 : end-4)];
            print(hf, '-djpeg', [path_name_psd, '.jpg'], '-noui');
            % Save to text file
            writetable( ...
                table( ...
                    tbl_data.FrequencyVector{1}, psd, SPL, ...
                    'VariableNames', {'freq, Hz', 'PSD, Pa^2/Hz', ['SPL' SPL_units_label]} ...
                    ), ...
                [path_name_psd '.csv'] ...
                )
        end
    end

    if false
    % Adaptive SLs (ASLs) adjust their order to the central frequency to
    % compensate the increasing wavelet bandwidth with increasing frequency.
    figure;
    subplot(1, 2, 1);
    imagesc(time, fois, aslt(xSignal, fs, fois, 3, srord, 0));
    colormap jet;

    %subplot(1, 1, 1);
    imagesc(time, fois, aslt(xSignal, fs, fois, 5, srord, 0));
    % The fractional ASLT (FASLT) provides sharp representations across
    % the entire frequency domain
    set(gca, 'ydir', 'normal');
    colormap jet;

    % Periodogram
    N = 128;
    T = 1;
    t = linspace(0,T,N);
    x = 12*sin(2*pi*10*t+pi/4)+5*cos(2*pi*40*t);
    figure,plot(x)
    dt = t(2)-t(1);
    f = 1/dt;
    X= fft(x);
    F = X(1:N/2+1);
    f = f*(0:N/2)/N;
    figure, plot(f,abs(F),'-');
    xlabel('Frequency');
    ylabel('|F(k)|');
    end
end

function [ar, fs, tim] = read_data(filename, coef)
    % Gets input array, its frequncy in Hz and recording time from text file
    %
    % filename: string, full path
    % coef: ar multiplier coefficient
    % Returns:
    % ar: 1D double, vector
    % fs: input rate, Hz
    % tim: datetime, time of data recording
    %
    % %%% File format %%%
    % Oscilloscope Data File
    % Experiment Time :   25-08-2022 21:25:32
    % Number of frames: 1639889
    % 
    % Module: E-502 (3T439094)
    % 
    % Number Of Channels : 1
    % Input Rate In kHz: 100.000000
    % Input Time In Sec: 16.398890
    % Decimation: 1
    % Data Format: Volts
    % Segments: 1
    % Data as Time Sequence:
    %     Ch  1  
    %  Канал 1   
    % 
    %      0.016
    %      0.012
    %      ...

    nl = '\r\n';
    fileID = fopen(filename, 'r', 'native', 'ASCII');
        fgetl(fileID);  % skip 1 line
        tim_str = fgetl(fileID);
        try
            str_arr = split(tim_str, ': ');
            tim_str = strtrim(str_arr{2});
            tim = datetime(tim_str, 'InputFormat', 'dd-MM-yyyy hh:mm:ss');
        catch
            tim = tim_str;
        end        
        headers = textscan(fileID, '%s %s', 9, 'Delimiter', [':' nl]);
        textscan(fileID, '%*s', 4, 'Delimiter', nl);
        ar = textscan(fileID, '%f', 'Delimiter', nl);
    fclose(fileID);
    ar = ar{1} * coef;
    fs = str2double(headers{1, 2}{strcmp(headers{1, 1}, ...
        'Input Rate In kHz')});
    fs = fs * 1E3;
end

function filterBankSA = get_spectrum_analyzer(fs)
    % Filter-bank
    
	global lim_x lim_y
    n_freq_bins = 1024*16;  % (for example 32000 is max for ZetLab)
    filterBankRBW = fs/n_freq_bins;
    try
        filterBankSA = spectrumAnalyzer( ...
            'SampleRate', fs, ...
            'RBWSource', 'property', ...  % means using RBW value:
            'RBW', filterBankRBW, ...     % Resolution bandwidth, Hz. (frequency span/RBW > 2)
            ... % 'AveragingMethod','exponential', ...
            ... % 'ForgettingFactor',0.001, ...
            'PlotAsTwoSidedSpectrum', false, ...
            'FrequencyScale', 'log', ...
            'SpectrumType', 'rms', ...  % 'power',  ... %
            'SpectrumUnits', 'Vrms', ...  % 'Watts', ... %
            ...  % YLabel='Power', ...
            'Title', sprintf('Filter bank (of size %d) Power Spectrum Estimate', n_freq_bins), ...
            'FrequencySpan', 'start-and-stop-frequencies', ...
            'StartFrequency', lim_x(1), ...
            'StopFrequency', fs/2, ...
            ...  % YLimits=lim_y, ...  % 'YLimits',[-150 50],
            'Position', [50 375 800 450] ...
        );
    catch ME1  % Use of DSP System Toolbox function error. This file in stack:
        i_this = find(strcmp({ME1.stack.file}, ME1.stack(end).file), 1, 'first');
        fprintf(2, ...
            [' Skipped high resolution spectrum calculation because of ' ...
            '<a href="matlab:opentoline(''%s'',%d)">error</a>: "%s"\n'], ...
            ME1.stack(i_this).file, ME1.stack(i_this).line(1),  ME1.message);
        filterBankSA = [];
    end
end

function poctave_localplot_10lg(ha, P, CF)
    % Draw 10*log10(P).
    % ha: axes to draw on
    %
    % Changes comparative to ~ same code of MATLAB poctave(): 
    % - units: Hz only, standard frequencies
    % - added ha

    global lim_y
    PdB = 10*log10(P);
    % Determine engineering units
    % Changed: 0 instead of CF(end) to get Hz
    [~, scaleFactor, unitsStr] = signal.internal.utilities.getFrequencyEngUnits(0);
    % Plotting
    centerFreqLabels = categorical(scaleFactor.*to_standard_freq(CF));
    pObject = bar(ha, centerFreqLabels, PdB);
    pObject(1).BaseValue = lim_y(1);
    ylim(ha, lim_y);
    % Resolve the number of x-ticks to maxNumTicks.
    maxNumTicks = 10;
    if numel(CF)>maxNumTicks
        spc = ceil(numel(CF) / maxNumTicks);
        % If GUI is a bar plot, then xTicks should be a categorical array.
        xticks(ha, centerFreqLabels(1:spc:end))
    end
    xlabel(ha, [getString(message('signal:poctave:Frequency')) ', ' unitsStr])
end

function standard_freq = to_standard_freq(freqs)
    % Convert exact 1/3 octave frequencies to nearest that of ГОСТ 12090-80
    round2half = round(freqs*2, 2, 'significant')/2;
    round2 = round(freqs, 2, 'significant');
    round2_better = (abs(freqs - round2) - abs(freqs - round2half))./freqs < -0.002;
    standard_freq = round2half;
    standard_freq(round2_better) = round2(round2_better);
end