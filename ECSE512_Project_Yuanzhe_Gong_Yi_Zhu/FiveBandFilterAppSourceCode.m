classdef FiveBandFilterAppSourceCode < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        FiveBandDRRSFilter              matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        ResetButton                     matlab.ui.control.Button
        LowPassFilterCutoffFrequecyEditField  matlab.ui.control.NumericEditField
        LowPassFilterCutoffFrequecyEditFieldLabel  matlab.ui.control.Label
        Passband1CutoffFrequecyEditField  matlab.ui.control.NumericEditField
        Passband1CutoffFrequecyLabel    matlab.ui.control.Label
        Passband2CutoffFrequecyEditField  matlab.ui.control.NumericEditField
        Passband2CutoffFrequecyLabel    matlab.ui.control.Label
        HighPassFilterCutoffFrequecyEditField  matlab.ui.control.NumericEditField
        HighPassFilterCutoffFrequecyEditFieldLabel  matlab.ui.control.Label
        LowPassFilterGainEditField      matlab.ui.control.NumericEditField
        LowPassFilterGainEditFieldLabel  matlab.ui.control.Label
        Bandpass1GainEditField          matlab.ui.control.NumericEditField
        Bandpass1GainLabel              matlab.ui.control.Label
        Bandpass2GainEditField          matlab.ui.control.NumericEditField
        Bandpass2GainLabel              matlab.ui.control.Label
        Bandpass3GainEditField          matlab.ui.control.NumericEditField
        Bandpass3GainEditFieldLabel     matlab.ui.control.Label
        HighPassFilterGainEditField     matlab.ui.control.NumericEditField
        HighPassFilterGainEditFieldLabel  matlab.ui.control.Label
        ThesamplingfrquencyofyourfileisEditField  matlab.ui.control.NumericEditField
        ThesamplingfrquencyofyourfileisLabel  matlab.ui.control.Label
        PleaseenterthefullnameofyouraudiofileEditField  matlab.ui.control.EditField
        PleaseenterthefullnameofyouraudiofileEditFieldLabel  matlab.ui.control.Label
        HzLabel                         matlab.ui.control.Label
        ImportButton                    matlab.ui.control.Button
        PleaseenterthefollowingparameterstocontrolyouraudioLabel  matlab.ui.control.Label
        ParameterSettingLabel           matlab.ui.control.Label
        ConfirmButton                   matlab.ui.control.Button
        HzLabel_2                       matlab.ui.control.Label
        HzLabel_3                       matlab.ui.control.Label
        HzLabel_4                       matlab.ui.control.Label
        HzLabel_5                       matlab.ui.control.Label
        dBLabel                         matlab.ui.control.Label
        dBLabel_2                       matlab.ui.control.Label
        dBLabel_3                       matlab.ui.control.Label
        dBLabel_4                       matlab.ui.control.Label
        dBLabel_5                       matlab.ui.control.Label
        RightPanel                      matlab.ui.container.Panel
        XaxislimitHzSliderLabel         matlab.ui.control.Label
        XaxislimitHzSlider              matlab.ui.control.Slider
        DualRecursiveRunningSumToneControlSystemLabel  matlab.ui.control.Label
        ChooseplotcontentDropDownLabel  matlab.ui.control.Label
        ChooseplotcontentDropDown       matlab.ui.control.DropDown
        ExportandPlotButton             matlab.ui.control.Button
        SavetheFigureButton             matlab.ui.control.Button
        YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel  matlab.ui.control.Label
        ReadytosavetheresultLampLabel   matlab.ui.control.Label
        ReadytosavetheresultLamp        matlab.ui.control.Lamp
        CustomizeaxislimitCheckBox      matlab.ui.control.CheckBox
        YaxislimitSliderLabel           matlab.ui.control.Label
        YaxislimitSlider                matlab.ui.control.Slider
        OutputaudiofileEditField        matlab.ui.control.EditField
        OutputaudiofileEditFieldLabel   matlab.ui.control.Label
        dBScaleCheckBox                 matlab.ui.control.CheckBox
        UIAxes                          matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    methods (Access = private)
        
        function mainfunction=inputparam(app, FLC,FB1,FB2,FHC,GL,GB1,GB2,GB3,GH,FS)
            syms z
            L1=app.solveLength(FLC,FS);
            B1=app.solveLength(FB1,FS);
            B11=app.solveLength(FB2,FS);
            H1=app.solveLength(FHC,FS);
            mainfunction=app.Equalizer(L1,B1,B11,H1,GL,GB1,GB2,GB3,GH);
        end
        
        function result1=Drrs(app, L1,L2)
            syms z
            result1=(1-z^(-L1))/(1-z^(-1))*(1-z^(-L2))/(1-z^(-1));
        end
        
        function y=getodd(app, x)
            y = 2*floor(x/2)+1;
        end
        
        function bandLength=solveLength(app,fc,fs)
            bandLength=app.getodd(1/(fc*2/fs));
        end
        
        function result2=Equalizer(app, L1,B1,B11,H1,GL,GB1,GB2,GB3,GH)
            syms z
            L2=app.getodd(L1/sqrt(2));
            H2=app.getodd(H1/sqrt(2));
            B2=app.getodd(B1/sqrt(2));
            B22=app.getodd(B11/sqrt(2));
            kH=1/(H1*H2);
            kL=1/(L1*L2);
            k1=1/(B1*B2);
            k2=1/(B11*B22);
            GL=10^(GL/20);
            GB1=10^(GB1/20);
            GB2=10^(GB2/20);
            GB3=10^(GB3/20);
            GH=10^(GH/20);
            
            result2=GH*z^(-[L1+L2+B1+B2+B11+B22]/2+1)*[z^(-[H1+H2]/2+1)-kH...
                *app.Drrs(H1,H2)]+[GB3]*[kH*z^(-[L1+L2+B1+B2+B11+B22]/2+1)...
                *app.Drrs(H1,H2)-k2*z^(-[H1+H2+B1+B2+L1+L2]/2+1)...
                *app.Drrs(B11,B22)]+[GB2]*[k2*z^(-[H1+H2+B1+B2+L1+L2]/2+1)...
                *app.Drrs(B11,B22)-k1*z^(-[H1+H2+L1+L2+B11+B22]/2+1)...
                *app.Drrs(B1,B2)]+(GB1)*[k1*z^(-[H1+H2+L1+L2+B11+B22]/2+1)...
                *app.Drrs(B1,B2)-kL*z^(-[H1+H2+B1+B2+B11+B22]/2+1)...
                *app.Drrs(L1,L2)]+GL*kL*z^(-[H1+H2+B1+B2+B11+B22]/2+1)...
                *app.Drrs(L1,L2);
         end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.FiveBandDRRSFilter.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {678, 678};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {382, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Button pushed function: ExportandPlotButton
        function ExportandPlotButtonPushed(app, event)
            cla(app.UIAxes);
            app.ResetButton.BackgroundColor = [0.00,0.45,0.74];
            app.ResetButton.FontColor = [1,1,1];
            
            % Get input file and display sampling frequency
            value = app.PleaseenterthefullnameofyouraudiofileEditField.Value;
            [x,Fs] = audioread(value);
            app.ThesamplingfrquencyofyourfileisEditField.Value = Fs;
            
            % Get input values
            flc = app.LowPassFilterCutoffFrequecyEditField.Value;
            fb1 = app.Passband1CutoffFrequecyEditField.Value;
            fb2 = app.Passband2CutoffFrequecyEditField.Value;
            fhc = app.HighPassFilterCutoffFrequecyEditField.Value;
            gl = app.LowPassFilterGainEditField.Value;
            gb1 = app.Bandpass1GainEditField.Value;
            gb2 = app.Bandpass2GainEditField.Value;
            gb3 = app.Bandpass3GainEditField.Value;
            gh = app.HighPassFilterGainEditField.Value;
            xlimit = app.XaxislimitHzSlider.Value;
            ylimit = app.YaxislimitSlider.Value;
            cus = app.CustomizeaxislimitCheckBox.Value;
            filename1 = app.OutputaudiofileEditField.Value;
            dB = app.dBScaleCheckBox.Value;
            syms z
            
            % Input parameter
            % inputparam(FLC,FB1,FB2,FHC,GL,GB1,GB2,GB3,GH,FS)
            system=app.inputparam(flc, fb1, fb2, fhc, gl, gb1, gb2, gb3,...
                gh, Fs); 
            
            % Designed filter based on pass band cutoff and target gain
            % results=simplify(system); 
            [fzn, fzd] = numden(system);
            fzc = sym2poly(fzn);
            fzb = sym2poly(fzd);
            
            % Get results after filtering
            y = filter(fzc,fzb,x);
            
            % Time domain impulse response convolution method
            %h=impz(fzc,fzb);
            %fil1=conv(x,h);
            %y=fil1(length(h):length(fil1));
            audiowrite(filename1,y,Fs);
            
            switch app.ChooseplotcontentDropDown.Value
                
                % Plot original signal before filtering
                % Plot in freqeuncy domain
                case 'Original Signal'
                    NFFT=length(x);
                    X=fft(x,NFFT);
                    F = ((0:1/NFFT:1-1/NFFT)*Fs).';
                    magnitudeX = abs(X);        % Magnitude of the FFT
                    %phaseX = unwrap(angle(X));  % Phase of the FFT
                    if (dB == 1)
                        if (cus == 1)
                            xlim(app.UIAxes,[0 xlimit])
                            ylim(app.UIAxes,[0 ylimit])
                            p = plot(app.UIAxes, F, 20*log10(magnitudeX), 'Color', [0.00,0.45,0.74]);
                            legend(app.UIAxes, p ,{'Original Signal'})
                        elseif (cus == 0)
                            xlim(app.UIAxes, [0 Fs/2])
                            ylim(app.UIAxes, 'auto')
                            p = plot(app.UIAxes, F, 20*log10(magnitudeX), 'Color', [0.00,0.45,0.74]);
                            legend(app.UIAxes, p ,{'Original Signal'})
                        end
                    elseif (dB == 0)
                        if (cus == 1)
                            xlim(app.UIAxes,[0 xlimit])
                            ylim(app.UIAxes,[0 ylimit])
                            p = plot(app.UIAxes, F, magnitudeX, 'Color', [0.00,0.45,0.74]);
                            legend(app.UIAxes, p ,{'Original Signal'})
                        elseif (cus == 0)
                            xlim(app.UIAxes, [0 Fs/2])
                            ylim(app.UIAxes, 'auto')
                            p = plot(app.UIAxes, F, magnitudeX, 'Color', [0.00,0.45,0.74]);
                            legend(app.UIAxes, p ,{'Original Signal'})
                        end
                    end
                    app.ReadytosavetheresultLamp.Color = [0 1 0];
                    
                % Plot original signal after filtering
                % Plot in freqeuncy domain    
                case 'Filtered Signal'
                    NFFT=length(y);
                    Y=fft(y,NFFT);
                    F = ((0:1/NFFT:1-1/NFFT)*Fs).';
                    magnitudeY = abs(Y);        % Magnitude of the FFT
                    %phaseY = unwrap(angle(Y));  % Phase of the FFT 
                    if (dB == 1)
                        %plot(app.UIAxes, F, magnitudeY)
                        % Customize axis limit
                        if (cus == 1)
                            xlim(app.UIAxes,[0 xlimit])
                            ylim(app.UIAxes,[0 ylimit])
                            p = plot(app.UIAxes, F, 20*log10(magnitudeY), 'Color', [0.85,0.325,0.098]);
                            legend(app.UIAxes, p ,{'Filtered Signal'})
                        elseif (cus == 0)
                            xlim(app.UIAxes, [0 Fs/2])
                            ylim(app.UIAxes, 'auto')
                            p = plot(app.UIAxes, F, 20*log10(magnitudeY), 'Color', [0.85,0.325,0.098]);
                            legend(app.UIAxes, p, {'Filtered Signal'})
                        end
                  elseif (dB == 0)
                        % Customize axis limit
                        if (cus == 1)
                            xlim(app.UIAxes,[0 xlimit])
                            ylim(app.UIAxes,[0 ylimit])
                            p = plot(app.UIAxes, F, magnitudeY, 'Color', [0.85,0.325,0.098]);
                            legend(app.UIAxes, p ,{'Filtered Signal'})
                        elseif (cus == 0)
                            xlim(app.UIAxes, [0 Fs/2])
                            ylim(app.UIAxes, 'auto')
                            p = plot(app.UIAxes, F, magnitudeY, 'Color', [0.85,0.325,0.098]);
                            legend(app.UIAxes, p, {'Filtered Signal'})
                        end
                    end
                    app.ReadytosavetheresultLamp.Color = [0 1 0];
                
                % Plot both signals
                % Plot in freqeuncy domain 
                case 'Both Signals'
                    % plot everything in frequency domain
                    NFFT=length(x);
                    X=fft(x,NFFT);
                    F = ((0:1/NFFT:1-1/NFFT)*Fs).';
                    magnitudeX = abs(X);
                    %phaseX = unwrap(angle(X));
                 
                    ax0=app.UIAxes;
                    hold(ax0, 'on')
                    
                    % Plot after filtering data
                    NFFT=length(y);
                    Y=fft(y,NFFT);
                    F = ((0:1/NFFT:1-1/NFFT)*Fs).';
                    magnitudeY = abs(Y);
                    %phaseY = unwrap(angle(Y)); 
                    
                    % Customize axis limit
                    if (dB == 0)
                        if (cus == 1)
                            xlim(app.UIAxes,[0 xlimit])
                            ylim(app.UIAxes,[0 ylimit])
                            p1=plot(app.UIAxes,F,magnitudeY, 'Color', [0.85,0.325,0.098]);
                            hold(ax0, 'on')
                            p2=plot(app.UIAxes, F, magnitudeX, 'Color', [0.00,0.45,0.74]);
                            hold(ax0, 'off')
                            legend(app.UIAxes, [p1 p2],{'Filtered Signal', 'Original Signal'})
                        elseif (cus == 0)
                            xlim(app.UIAxes, [0 Fs/2])
                            ylim(app.UIAxes, 'auto')
                            p1=plot(app.UIAxes,F,magnitudeY, 'Color', [0.85,0.325,0.098]);
                            hold(ax0, 'on')
                            p2=plot(app.UIAxes, F, magnitudeX, 'Color', [0.00,0.45,0.74]);
                            hold(ax0, 'off')
                            legend(app.UIAxes, [p1 p2],{'Filtered Signal', 'Original Signal'})
                        end
                    elseif (dB == 1)
                        if (cus == 1)
                            xlim(app.UIAxes,[0 xlimit])
                            ylim(app.UIAxes,[0 ylimit])
                            p1=plot(app.UIAxes,F,20*log10(magnitudeY), 'Color', [0.85,0.325,0.098]);
                            hold(ax0, 'on')
                            p2=plot(app.UIAxes, F, 20*log10(magnitudeX), 'Color', [0.00,0.45,0.74]);
                            hold(ax0, 'off')
                            legend(app.UIAxes, [p1 p2],{'Filtered Signal', 'Original Signal'})
                        elseif (cus == 0)
                            xlim(app.UIAxes, [0 Fs/2])
                            ylim(app.UIAxes, 'auto')
                            p1=plot(app.UIAxes,F, 20*log10(magnitudeY), 'Color', [0.85,0.325,0.098]);
                            hold(ax0, 'on')
                            p2=plot(app.UIAxes, F, 20*log10(magnitudeX), 'Color', [0.00,0.45,0.74]);
                            hold(ax0, 'off')
                            legend(app.UIAxes, [p1 p2],{'Filtered Signal', 'Original Signal'})
                        end
                        
                    end
                    
                    app.ReadytosavetheresultLamp.Color = [0 1 0];
            end 
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            % Reset all figures and parameters
            cla(app.UIAxes);
            app.ReadytosavetheresultLamp.Color = [0.65,0.65,0.65];
            app.ExportandPlotButton.BackgroundColor=[0.65,0.65,0.65];
            app.ImportButton.BackgroundColor=[0.00,0.45,0.74];
            app.XaxislimitHzSlider.Value= 1000;
            app.YaxislimitSlider.Value=100;
            app.ThesamplingfrquencyofyourfileisEditField.Value = 0;
            app.LowPassFilterCutoffFrequecyEditField.Value = 1;
            app.Passband1CutoffFrequecyEditField.Value = 1;
            app.Passband2CutoffFrequecyEditField.Value = 1;
            app.HighPassFilterCutoffFrequecyEditField.Value = 1;
            app.LowPassFilterCutoffFrequecyEditField.Value = 1;
            app.LowPassFilterGainEditField.Value = 0;
            app.Bandpass1GainEditField.Value = 0;
            app.Bandpass2GainEditField.Value = 0;
            app.Bandpass3GainEditField.Value = 0;
            app.HighPassFilterGainEditField.Value = 0;
            app.ResetButton.BackgroundColor = [0.96,0.96,0.96];
            app.ResetButton.FontColor=[0,0,0];
            app.CustomizeaxislimitCheckBox.Value = 0;
            app.dBScaleCheckBox.Value = 0;
            app.ChooseplotcontentDropDown.Value = 'Original Signal';
            
        end

        % Value changed function: 
        % LowPassFilterCutoffFrequecyEditField
        function LowPassFilterCutoffFrequecyEditFieldValueChanged(app, event)
            value = app.LowPassFilterCutoffFrequecyEditField.Value;
            
        end

        % Value changed function: Passband1CutoffFrequecyEditField
        function Passband1CutoffFrequecyEditFieldValueChanged(app, event)
            value = app.Passband1CutoffFrequecyEditField.Value;
            
        end

        % Value changed function: Passband2CutoffFrequecyEditField
        function Passband2CutoffFrequecyEditFieldValueChanged(app, event)
            value = app.Passband2CutoffFrequecyEditField.Value;
            
        end

        % Value changed function: 
        % HighPassFilterCutoffFrequecyEditField
        function HighPassFilterCutoffFrequecyEditFieldValueChanged(app, event)
            value = app.HighPassFilterCutoffFrequecyEditField.Value;
            
        end

        % Value changed function: LowPassFilterGainEditField
        function LowPassFilterGainEditFieldValueChanged(app, event)
            value = app.LowPassFilterGainEditField.Value;
            
        end

        % Value changed function: Bandpass1GainEditField
        function Bandpass1GainEditFieldValueChanged(app, event)
            value = app.Bandpass1GainEditField.Value;
            
        end

        % Value changed function: Bandpass2GainEditField
        function Bandpass2GainEditFieldValueChanged(app, event)
            value = app.Bandpass2GainEditField.Value;
            
        end

        % Value changed function: Bandpass3GainEditField
        function Bandpass3GainEditFieldValueChanged(app, event)
            value = app.Bandpass3GainEditField.Value;
            
        end

        % Value changed function: HighPassFilterGainEditField
        function HighPassFilterGainEditFieldValueChanged(app, event)
            value = app.HighPassFilterGainEditField.Value;
            
        end

        % Button pushed function: ImportButton
        function ImportButtonPushed(app, event)
            value = app.PleaseenterthefullnameofyouraudiofileEditField.Value;
            [x,Fs] = audioread(value);
            app.ThesamplingfrquencyofyourfileisEditField.Value = Fs;
            app.ConfirmButton.BackgroundColor=[0.00,0.45,0.74];
            app.ImportButton.BackgroundColor=[0.65,0.65,0.65];
        end

        % Value changed function: ChooseplotcontentDropDown
        function ChooseplotcontentDropDownValueChanged(app, event)
            value = app.ChooseplotcontentDropDown.Value;
        end

        % Button pushed function: SavetheFigureButton
        function SavetheFigureButtonPushed(app, event)
            % Create a temporary figure with axes.
            figure;
            % Create new axis
            newAxes = axes;
            
            % Copy all objects from UIAxes to new axis
            newAxes.Title= app.UIAxes.Title;
            newAxes.XLabel = app.UIAxes.XLabel;
            newAxes.YLabel = app.UIAxes.YLabel;
            copyobj(app.UIAxes.Children, newAxes)
            
            % Keep titles and labels on uiaxes
            xlabel(app.UIAxes, newAxes.XLabel.String);
            ylabel(app.UIAxes, newAxes.YLabel.String);
            title(app.UIAxes, newAxes.Title.String);
            
            fig2Save = ancestor(newAxes, 'figure');
            fig2Save.Visible='off';
           
            % Save as window
            filter = {'*.jpg';'*.png'};
            [filename, pathname] = uiputfile(filter);
            newfilename = fullfile(pathname, filename);
            saveas(fig2Save, newfilename);
        end

        % Button pushed function: ConfirmButton
        function ConfirmButtonPushed(app, event)
            app.ExportandPlotButton.BackgroundColor=[0.00,0.45,0.74];
            app.ConfirmButton.BackgroundColor=[0.65,0.65,0.65];
        end

        % Value changed function: XaxislimitHzSlider
        function XaxislimitHzSliderValueChanged(app, event)
            value = app.XaxislimitHzSlider.Value;
            
        end

        % Value changed function: CustomizeaxislimitCheckBox
        function CustomizeaxislimitCheckBoxValueChanged(app, event)
            value = app.CustomizeaxislimitCheckBox.Value;
            
        end

        % Value changed function: YaxislimitSlider
        function YaxislimitSliderValueChanged(app, event)
            value = app.YaxislimitSlider.Value;
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create FiveBandDRRSFilter and hide until all components are created
            app.FiveBandDRRSFilter = uifigure('Visible', 'off');
            app.FiveBandDRRSFilter.AutoResizeChildren = 'off';
            app.FiveBandDRRSFilter.Position = [100 100 1021 678];
            app.FiveBandDRRSFilter.Name = 'FiveBandDRRSFilter';
            app.FiveBandDRRSFilter.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.FiveBandDRRSFilter);
            app.GridLayout.ColumnWidth = {382, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Scrollable = 'on';

            % Create ResetButton
            app.ResetButton = uibutton(app.LeftPanel, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [69 34 100 22];
            app.ResetButton.Text = 'Reset';

            % Create LowPassFilterCutoffFrequecyEditField
            app.LowPassFilterCutoffFrequecyEditField = uieditfield(app.LeftPanel, 'numeric');
            app.LowPassFilterCutoffFrequecyEditField.Limits = [1 10000000];
            app.LowPassFilterCutoffFrequecyEditField.ValueDisplayFormat = '%8.f';
            app.LowPassFilterCutoffFrequecyEditField.ValueChangedFcn = createCallbackFcn(app, @LowPassFilterCutoffFrequecyEditFieldValueChanged, true);
            app.LowPassFilterCutoffFrequecyEditField.Position = [229 407 100 22];
            app.LowPassFilterCutoffFrequecyEditField.Value = 300;

            % Create LowPassFilterCutoffFrequecyEditFieldLabel
            app.LowPassFilterCutoffFrequecyEditFieldLabel = uilabel(app.LeftPanel);
            app.LowPassFilterCutoffFrequecyEditFieldLabel.HorizontalAlignment = 'right';
            app.LowPassFilterCutoffFrequecyEditFieldLabel.Position = [34 407 180 22];
            app.LowPassFilterCutoffFrequecyEditFieldLabel.Text = 'Low Pass Filter Cut-off Frequecy';

            % Create Passband1CutoffFrequecyEditField
            app.Passband1CutoffFrequecyEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Passband1CutoffFrequecyEditField.Limits = [1 10000000];
            app.Passband1CutoffFrequecyEditField.ValueChangedFcn = createCallbackFcn(app, @Passband1CutoffFrequecyEditFieldValueChanged, true);
            app.Passband1CutoffFrequecyEditField.Position = [229 367 100 22];
            app.Passband1CutoffFrequecyEditField.Value = 900;

            % Create Passband1CutoffFrequecyLabel
            app.Passband1CutoffFrequecyLabel = uilabel(app.LeftPanel);
            app.Passband1CutoffFrequecyLabel.HorizontalAlignment = 'right';
            app.Passband1CutoffFrequecyLabel.Position = [46 367 168 22];
            app.Passband1CutoffFrequecyLabel.Text = 'Passband #1 Cut-off Frequecy';

            % Create Passband2CutoffFrequecyEditField
            app.Passband2CutoffFrequecyEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Passband2CutoffFrequecyEditField.Limits = [1 10000000];
            app.Passband2CutoffFrequecyEditField.ValueDisplayFormat = '%.0f';
            app.Passband2CutoffFrequecyEditField.ValueChangedFcn = createCallbackFcn(app, @Passband2CutoffFrequecyEditFieldValueChanged, true);
            app.Passband2CutoffFrequecyEditField.Position = [229 327 100 22];
            app.Passband2CutoffFrequecyEditField.Value = 1500;

            % Create Passband2CutoffFrequecyLabel
            app.Passband2CutoffFrequecyLabel = uilabel(app.LeftPanel);
            app.Passband2CutoffFrequecyLabel.HorizontalAlignment = 'right';
            app.Passband2CutoffFrequecyLabel.Position = [46 327 168 22];
            app.Passband2CutoffFrequecyLabel.Text = 'Passband #2 Cut-off Frequecy';

            % Create HighPassFilterCutoffFrequecyEditField
            app.HighPassFilterCutoffFrequecyEditField = uieditfield(app.LeftPanel, 'numeric');
            app.HighPassFilterCutoffFrequecyEditField.Limits = [1 10000000];
            app.HighPassFilterCutoffFrequecyEditField.ValueDisplayFormat = '%.0f';
            app.HighPassFilterCutoffFrequecyEditField.ValueChangedFcn = createCallbackFcn(app, @HighPassFilterCutoffFrequecyEditFieldValueChanged, true);
            app.HighPassFilterCutoffFrequecyEditField.Position = [229 287 100 22];
            app.HighPassFilterCutoffFrequecyEditField.Value = 2500;

            % Create HighPassFilterCutoffFrequecyEditFieldLabel
            app.HighPassFilterCutoffFrequecyEditFieldLabel = uilabel(app.LeftPanel);
            app.HighPassFilterCutoffFrequecyEditFieldLabel.HorizontalAlignment = 'right';
            app.HighPassFilterCutoffFrequecyEditFieldLabel.Position = [31 287 183 22];
            app.HighPassFilterCutoffFrequecyEditFieldLabel.Text = 'High Pass Filter Cut-off Frequecy';

            % Create LowPassFilterGainEditField
            app.LowPassFilterGainEditField = uieditfield(app.LeftPanel, 'numeric');
            app.LowPassFilterGainEditField.Limits = [-10000000 10000000];
            app.LowPassFilterGainEditField.ValueDisplayFormat = '%.0f';
            app.LowPassFilterGainEditField.ValueChangedFcn = createCallbackFcn(app, @LowPassFilterGainEditFieldValueChanged, true);
            app.LowPassFilterGainEditField.Position = [229 247 100 22];

            % Create LowPassFilterGainEditFieldLabel
            app.LowPassFilterGainEditFieldLabel = uilabel(app.LeftPanel);
            app.LowPassFilterGainEditFieldLabel.HorizontalAlignment = 'right';
            app.LowPassFilterGainEditFieldLabel.Position = [100 247 114 22];
            app.LowPassFilterGainEditFieldLabel.Text = 'Low Pass Filter Gain';

            % Create Bandpass1GainEditField
            app.Bandpass1GainEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Bandpass1GainEditField.Limits = [-10000000 10000000];
            app.Bandpass1GainEditField.ValueDisplayFormat = '%.0f';
            app.Bandpass1GainEditField.ValueChangedFcn = createCallbackFcn(app, @Bandpass1GainEditFieldValueChanged, true);
            app.Bandpass1GainEditField.Position = [229 207 100 22];
            app.Bandpass1GainEditField.Value = 10;

            % Create Bandpass1GainLabel
            app.Bandpass1GainLabel = uilabel(app.LeftPanel);
            app.Bandpass1GainLabel.HorizontalAlignment = 'right';
            app.Bandpass1GainLabel.Position = [100 207 114 22];
            app.Bandpass1GainLabel.Text = 'Bandpass #1 Gain';

            % Create Bandpass2GainEditField
            app.Bandpass2GainEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Bandpass2GainEditField.Limits = [-10000000 10000000];
            app.Bandpass2GainEditField.ValueDisplayFormat = '%.0f';
            app.Bandpass2GainEditField.ValueChangedFcn = createCallbackFcn(app, @Bandpass2GainEditFieldValueChanged, true);
            app.Bandpass2GainEditField.Position = [229 167 100 22];

            % Create Bandpass2GainLabel
            app.Bandpass2GainLabel = uilabel(app.LeftPanel);
            app.Bandpass2GainLabel.HorizontalAlignment = 'right';
            app.Bandpass2GainLabel.Position = [110 167 104 22];
            app.Bandpass2GainLabel.Text = 'Bandpass #2 Gain';

            % Create Bandpass3GainEditField
            app.Bandpass3GainEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Bandpass3GainEditField.Limits = [-10000000 10000000];
            app.Bandpass3GainEditField.ValueDisplayFormat = '%.0f';
            app.Bandpass3GainEditField.ValueChangedFcn = createCallbackFcn(app, @Bandpass3GainEditFieldValueChanged, true);
            app.Bandpass3GainEditField.Position = [229 127 100 22];

            % Create Bandpass3GainEditFieldLabel
            app.Bandpass3GainEditFieldLabel = uilabel(app.LeftPanel);
            app.Bandpass3GainEditFieldLabel.HorizontalAlignment = 'right';
            app.Bandpass3GainEditFieldLabel.Position = [110 127 104 22];
            app.Bandpass3GainEditFieldLabel.Text = 'Bandpass #3 Gain';

            % Create HighPassFilterGainEditField
            app.HighPassFilterGainEditField = uieditfield(app.LeftPanel, 'numeric');
            app.HighPassFilterGainEditField.Limits = [-10000000 10000000];
            app.HighPassFilterGainEditField.ValueDisplayFormat = '%.0f';
            app.HighPassFilterGainEditField.ValueChangedFcn = createCallbackFcn(app, @HighPassFilterGainEditFieldValueChanged, true);
            app.HighPassFilterGainEditField.Position = [229 87 100 22];

            % Create HighPassFilterGainEditFieldLabel
            app.HighPassFilterGainEditFieldLabel = uilabel(app.LeftPanel);
            app.HighPassFilterGainEditFieldLabel.HorizontalAlignment = 'right';
            app.HighPassFilterGainEditFieldLabel.Position = [95 87 119 22];
            app.HighPassFilterGainEditFieldLabel.Text = 'High Pass Filter Gain';

            % Create ThesamplingfrquencyofyourfileisEditField
            app.ThesamplingfrquencyofyourfileisEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ThesamplingfrquencyofyourfileisEditField.Limits = [0 10000000];
            app.ThesamplingfrquencyofyourfileisEditField.ValueDisplayFormat = '%.0f';
            app.ThesamplingfrquencyofyourfileisEditField.Position = [229 492 100 22];

            % Create ThesamplingfrquencyofyourfileisLabel
            app.ThesamplingfrquencyofyourfileisLabel = uilabel(app.LeftPanel);
            app.ThesamplingfrquencyofyourfileisLabel.Position = [20 492 197 22];
            app.ThesamplingfrquencyofyourfileisLabel.Text = 'The sampling frquency of your file is';

            % Create PleaseenterthefullnameofyouraudiofileEditField
            app.PleaseenterthefullnameofyouraudiofileEditField = uieditfield(app.LeftPanel, 'text');
            app.PleaseenterthefullnameofyouraudiofileEditField.HorizontalAlignment = 'right';
            app.PleaseenterthefullnameofyouraudiofileEditField.Position = [253 569 100 22];
            app.PleaseenterthefullnameofyouraudiofileEditField.Value = 'sample.wav';

            % Create PleaseenterthefullnameofyouraudiofileEditFieldLabel
            app.PleaseenterthefullnameofyouraudiofileEditFieldLabel = uilabel(app.LeftPanel);
            app.PleaseenterthefullnameofyouraudiofileEditFieldLabel.HorizontalAlignment = 'center';
            app.PleaseenterthefullnameofyouraudiofileEditFieldLabel.Position = [4 569 250 22];
            app.PleaseenterthefullnameofyouraudiofileEditFieldLabel.Text = 'Please enter the full name of your audio file';

            % Create HzLabel
            app.HzLabel = uilabel(app.LeftPanel);
            app.HzLabel.HorizontalAlignment = 'right';
            app.HzLabel.Position = [328 492 25 22];
            app.HzLabel.Text = 'Hz';

            % Create ImportButton
            app.ImportButton = uibutton(app.LeftPanel, 'push');
            app.ImportButton.ButtonPushedFcn = createCallbackFcn(app, @ImportButtonPushed, true);
            app.ImportButton.BackgroundColor = [0 0.451 0.7412];
            app.ImportButton.FontColor = [1 1 1];
            app.ImportButton.Position = [142 533 100 22];
            app.ImportButton.Text = 'Import';

            % Create PleaseenterthefollowingparameterstocontrolyouraudioLabel
            app.PleaseenterthefollowingparameterstocontrolyouraudioLabel = uilabel(app.LeftPanel);
            app.PleaseenterthefollowingparameterstocontrolyouraudioLabel.HorizontalAlignment = 'center';
            app.PleaseenterthefollowingparameterstocontrolyouraudioLabel.Position = [6 447 356 22];
            app.PleaseenterthefollowingparameterstocontrolyouraudioLabel.Text = 'Please enter the following parameters to control your audio';

            % Create ParameterSettingLabel
            app.ParameterSettingLabel = uilabel(app.LeftPanel);
            app.ParameterSettingLabel.HorizontalAlignment = 'center';
            app.ParameterSettingLabel.FontName = 'Copperplate Gothic Bold';
            app.ParameterSettingLabel.FontSize = 20;
            app.ParameterSettingLabel.Position = [105 609 167 24];
            app.ParameterSettingLabel.Text = 'Parameter Setting';

            % Create ConfirmButton
            app.ConfirmButton = uibutton(app.LeftPanel, 'push');
            app.ConfirmButton.ButtonPushedFcn = createCallbackFcn(app, @ConfirmButtonPushed, true);
            app.ConfirmButton.BackgroundColor = [0.651 0.651 0.651];
            app.ConfirmButton.FontColor = [1 1 1];
            app.ConfirmButton.Position = [229 34 100 22];
            app.ConfirmButton.Text = 'Confirm';

            % Create HzLabel_2
            app.HzLabel_2 = uilabel(app.LeftPanel);
            app.HzLabel_2.HorizontalAlignment = 'right';
            app.HzLabel_2.Position = [328 407 25 22];
            app.HzLabel_2.Text = 'Hz';

            % Create HzLabel_3
            app.HzLabel_3 = uilabel(app.LeftPanel);
            app.HzLabel_3.HorizontalAlignment = 'right';
            app.HzLabel_3.Position = [328 367 25 22];
            app.HzLabel_3.Text = 'Hz';

            % Create HzLabel_4
            app.HzLabel_4 = uilabel(app.LeftPanel);
            app.HzLabel_4.HorizontalAlignment = 'right';
            app.HzLabel_4.Position = [328 327 25 22];
            app.HzLabel_4.Text = 'Hz';

            % Create HzLabel_5
            app.HzLabel_5 = uilabel(app.LeftPanel);
            app.HzLabel_5.HorizontalAlignment = 'right';
            app.HzLabel_5.Position = [328 287 25 22];
            app.HzLabel_5.Text = 'Hz';

            % Create dBLabel
            app.dBLabel = uilabel(app.LeftPanel);
            app.dBLabel.HorizontalAlignment = 'right';
            app.dBLabel.Position = [328 247 25 22];
            app.dBLabel.Text = 'dB';

            % Create dBLabel_2
            app.dBLabel_2 = uilabel(app.LeftPanel);
            app.dBLabel_2.HorizontalAlignment = 'right';
            app.dBLabel_2.Position = [328 207 25 22];
            app.dBLabel_2.Text = 'dB';

            % Create dBLabel_3
            app.dBLabel_3 = uilabel(app.LeftPanel);
            app.dBLabel_3.HorizontalAlignment = 'right';
            app.dBLabel_3.Position = [328 167 25 22];
            app.dBLabel_3.Text = 'dB';

            % Create dBLabel_4
            app.dBLabel_4 = uilabel(app.LeftPanel);
            app.dBLabel_4.HorizontalAlignment = 'right';
            app.dBLabel_4.Position = [328 127 25 22];
            app.dBLabel_4.Text = 'dB';

            % Create dBLabel_5
            app.dBLabel_5 = uilabel(app.LeftPanel);
            app.dBLabel_5.HorizontalAlignment = 'right';
            app.dBLabel_5.Position = [328 87 25 22];
            app.dBLabel_5.Text = 'dB';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Scrollable = 'on';

            % Create XaxislimitHzSliderLabel
            app.XaxislimitHzSliderLabel = uilabel(app.RightPanel);
            app.XaxislimitHzSliderLabel.HorizontalAlignment = 'center';
            app.XaxislimitHzSliderLabel.Position = [197 200 89 41];
            app.XaxislimitHzSliderLabel.Text = 'X-axis limit (Hz)';

            % Create XaxislimitHzSlider
            app.XaxislimitHzSlider = uislider(app.RightPanel);
            app.XaxislimitHzSlider.Limits = [1000 10000];
            app.XaxislimitHzSlider.ValueChangedFcn = createCallbackFcn(app, @XaxislimitHzSliderValueChanged, true);
            app.XaxislimitHzSlider.Position = [299 224 309 3];
            app.XaxislimitHzSlider.Value = 1000;

            % Create DualRecursiveRunningSumToneControlSystemLabel
            app.DualRecursiveRunningSumToneControlSystemLabel = uilabel(app.RightPanel);
            app.DualRecursiveRunningSumToneControlSystemLabel.HorizontalAlignment = 'center';
            app.DualRecursiveRunningSumToneControlSystemLabel.FontName = 'Copperplate Gothic Bold';
            app.DualRecursiveRunningSumToneControlSystemLabel.FontSize = 26;
            app.DualRecursiveRunningSumToneControlSystemLabel.FontWeight = 'bold';
            app.DualRecursiveRunningSumToneControlSystemLabel.Position = [6 628 645 33];
            app.DualRecursiveRunningSumToneControlSystemLabel.Text = 'Dual Recursive Running Sum Tone Control System';

            % Create ChooseplotcontentDropDownLabel
            app.ChooseplotcontentDropDownLabel = uilabel(app.RightPanel);
            app.ChooseplotcontentDropDownLabel.Position = [45 116 127 22];
            app.ChooseplotcontentDropDownLabel.Text = 'Choose plot content';

            % Create ChooseplotcontentDropDown
            app.ChooseplotcontentDropDown = uidropdown(app.RightPanel);
            app.ChooseplotcontentDropDown.Items = {'Original Signal', 'Filtered Signal', 'Both Signals'};
            app.ChooseplotcontentDropDown.ValueChangedFcn = createCallbackFcn(app, @ChooseplotcontentDropDownValueChanged, true);
            app.ChooseplotcontentDropDown.Position = [160 116 118 22];
            app.ChooseplotcontentDropDown.Value = 'Original Signal';

            % Create ExportandPlotButton
            app.ExportandPlotButton = uibutton(app.RightPanel, 'push');
            app.ExportandPlotButton.ButtonPushedFcn = createCallbackFcn(app, @ExportandPlotButtonPushed, true);
            app.ExportandPlotButton.BackgroundColor = [0.651 0.651 0.651];
            app.ExportandPlotButton.FontColor = [1 1 1];
            app.ExportandPlotButton.Position = [139 34 100 22];
            app.ExportandPlotButton.Text = 'Export and Plot';

            % Create SavetheFigureButton
            app.SavetheFigureButton = uibutton(app.RightPanel, 'push');
            app.SavetheFigureButton.ButtonPushedFcn = createCallbackFcn(app, @SavetheFigureButtonPushed, true);
            app.SavetheFigureButton.Position = [414 34 100 22];
            app.SavetheFigureButton.Text = 'Save the Figure';

            % Create YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel
            app.YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel = uilabel(app.RightPanel);
            app.YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel.HorizontalAlignment = 'center';
            app.YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel.FontName = 'Copperplate Gothic Bold';
            app.YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel.FontColor = [0.651 0.651 0.651];
            app.YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel.Position = [139 1 375 22];
            app.YuanzheGongYiZhuECEDepartmentMcGillUniversityLabel.Text = '© 2020 Yuanzhe Gong & Yi Zhu, ECE Department, McGill University';

            % Create ReadytosavetheresultLampLabel
            app.ReadytosavetheresultLampLabel = uilabel(app.RightPanel);
            app.ReadytosavetheresultLampLabel.Position = [387 95 138 22];
            app.ReadytosavetheresultLampLabel.Text = 'Ready to save the result';

            % Create ReadytosavetheresultLamp
            app.ReadytosavetheresultLamp = uilamp(app.RightPanel);
            app.ReadytosavetheresultLamp.Position = [524 96 20 20];
            app.ReadytosavetheresultLamp.Color = [0.651 0.651 0.651];

            % Create CustomizeaxislimitCheckBox
            app.CustomizeaxislimitCheckBox = uicheckbox(app.RightPanel);
            app.CustomizeaxislimitCheckBox.ValueChangedFcn = createCallbackFcn(app, @CustomizeaxislimitCheckBoxValueChanged, true);
            app.CustomizeaxislimitCheckBox.Text = 'Customize axis limit';
            app.CustomizeaxislimitCheckBox.Position = [44 162 128 22];

            % Create YaxislimitSliderLabel
            app.YaxislimitSliderLabel = uilabel(app.RightPanel);
            app.YaxislimitSliderLabel.HorizontalAlignment = 'center';
            app.YaxislimitSliderLabel.Position = [196 153 81 41];
            app.YaxislimitSliderLabel.Text = 'Y-axis limit';

            % Create YaxislimitSlider
            app.YaxislimitSlider = uislider(app.RightPanel);
            app.YaxislimitSlider.Limits = [100 10000];
            app.YaxislimitSlider.ValueChangedFcn = createCallbackFcn(app, @YaxislimitSliderValueChanged, true);
            app.YaxislimitSlider.Position = [299 177 309 3];
            app.YaxislimitSlider.Value = 100;

            % Create OutputaudiofileEditField
            app.OutputaudiofileEditField = uieditfield(app.RightPanel, 'text');
            app.OutputaudiofileEditField.HorizontalAlignment = 'right';
            app.OutputaudiofileEditField.Position = [158 76 119 22];
            app.OutputaudiofileEditField.Value = 'filteredsignal.wav';

            % Create OutputaudiofileEditFieldLabel
            app.OutputaudiofileEditFieldLabel = uilabel(app.RightPanel);
            app.OutputaudiofileEditFieldLabel.HorizontalAlignment = 'center';
            app.OutputaudiofileEditFieldLabel.Position = [44 76 95 22];
            app.OutputaudiofileEditFieldLabel.Text = 'Output audio file:';

            % Create dBScaleCheckBox
            app.dBScaleCheckBox = uicheckbox(app.RightPanel);
            app.dBScaleCheckBox.Text = 'dB Scale';
            app.dBScaleCheckBox.Position = [45 207 70 22];

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, 'Frequency Domain Signal Plot')
            xlabel(app.UIAxes, 'Frequency (Hz)')
            ylabel(app.UIAxes, 'Magnitude')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [25 250 600 364];

            % Show the figure after all components are created
            app.FiveBandDRRSFilter.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FiveBandFilterAppSourceCode

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.FiveBandDRRSFilter)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.FiveBandDRRSFilter)
        end
    end
end