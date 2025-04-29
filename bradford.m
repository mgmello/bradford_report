clear all; close all; clc
format long

% setting concentration values for albumine in standard curve

albu_concentration=0.2;

albu_volumes=[0 10 20 30 40 50 60 70 80];

final_volume=1.1;

albu_standard=albu_concentration*albu_volumes*1e-6/1e-3/final_volume;

max_albu=max(albu_standard);

% data processing

white595_1=[0.000 0.111 0.162 0.193 0.258 0.263 0.334 0.376 0.407];
white595_2=[0.000 0.101 0.149 0.185 0.251 0.291 0.339 0.352 0.383];
data.white595_mean=mean([white595_1; white595_2]);

di595_1=[0.342 0.262 0.068];
di595_2=[0.345 0.266 0.062];
sample.di595_mean=mean([di595_1; di595_2]);

dc595_1=[0.332 0.291 0.278];
dc595_2=[0.344 0.272 0.223];
sample.dc595_mean=mean([dc595_1; dc595_2]);

white540_1=[0.000 0.025 0.048 0.044 0.076 0.054 0.084 0.097 0.093];
white540_2=[0.000 0.017 0.034 0.039 0.078 0.070 0.091 0.083 0.080];
data.white540_mean=mean([white540_1; white540_2]);

di540_1=[0.000 0.014 0.014];
di540_2=[0.000 0.016 0.010];
sample.di540_mean=mean([di540_1; di540_2]);

dc540_1=[0.000 0.010 0.002];
dc540_2=[0.000 0.000 0.000];
sample.dc540_mean=mean([dc540_1; dc540_2]);

% best R2 and save linear coef

x_degree=3;
r2_dict=containers.Map();
r2_540_values_x=[];
r2_595_values_x=[];

for z=1:2
    if z==1
        wavelength='540';
    else
        wavelength='595';
    end
    for y=1:x_degree
        current_data_test=['white' wavelength '_mean'];
        current_group_test=data.(current_data_test);
        current_coef_r2_test=polyfit(albu_standard, current_group_test, y);
        y_pol_r2_x=polyval(current_coef_r2_test, albu_standard);
        sqs_r2_x=sum((current_group_test-y_pol_r2_x).^2);
        stq_r2_x=sum((current_group_test-mean(current_group_test)).^2);
        r2_test_x=1-(sqs_r2_x/stq_r2_x);
        if z==1
            r2_540_values_x(end+1)=r2_test_x;
            keyname=sprintf('540 degree %d', y);
            r2_dict(keyname)=r2_test_x;
        else
            r2_595_values_x(end+1)=r2_test_x;
            keyname=sprintf('595 degree %d', y);
            r2_dict(keyname)=r2_test_x;
        end
        if y==1
            if z==1
                lin.linear_coef_540=current_coef_r2_test;
            elseif z==2
                lin.linear_coef_595=current_coef_r2_test;
            end
        end
    end
end
[r2_540_max, idx540]=max(r2_540_values_x);
[r2_595_max, idx595]=max(r2_595_values_x);
for a=1:2
    if a==1
        keyhead='540 degree %d';
        r2_max=r2_540_max;
        best_degree=sprintf('degree %d', idx540);
        idx=idx540;
        wavelength='540';
    else
        keyhead='595 degree %d';
        r2_max=r2_595_max;
        best_degree=sprintf('degree %d', idx595);
        idx=idx595;
        wavelength='595';
    end
    fprintf('Best %s degree is %s and it R² value is %.20f\n', wavelength, best_degree, r2_max)
    for u=1:x_degree
        keyname=sprintf(keyhead, u);
        if u~=idx
            fprintf('For %s, R² value is %.20f\n', keyname, r2_dict(keyname))
        end
    end
end

% plotting bradford standard graph for 540nm and 595nm
value.coef_540=polyfit(albu_standard, data.white540_mean, idx540);
value.coef_595=polyfit(albu_standard, data.white595_mean, idx595);
x_fit=linspace(min(albu_standard), max(albu_standard), 100);
% R2 --------------------------------------------------------
y_pol_540=polyval(value.coef_540, albu_standard);
y_pol_595=polyval(value.coef_595, albu_standard);
sqs_540=sum((data.white540_mean-y_pol_540).^2);
sqs_595=sum((data.white595_mean-y_pol_595).^2);
stq_540=sum((data.white540_mean-mean(data.white540_mean)).^2);
stq_595=sum((data.white595_mean-mean(data.white595_mean)).^2);
R2_540=1-(sqs_540/stq_540);
R2_595=1-(sqs_595/stq_595);
%-------------------------------------------------------------
for m=1:2
    if m==1
        wavelength='540';
        idx_value=idx540;
    else
        wavelength='595';
        idx_value=idx595;
    end
    current_coef=['coef_' wavelength];
    current_poly=value.(current_coef);
    current_white=['white' wavelength '_mean'];
    current_data=data.(current_white);
    figure
    plot(albu_standard, current_data, 'o', 'LineWidth', 2)
    grid on
    xlabel('mg/mL', 'FontSize', 12)
    ylabel('abs', 'FontSize', 12)
    max_white_mean=max(current_data);
    ylim([0 max_white_mean+0.01])
    titulo=['Bradford standard curve - ' wavelength 'nm'];
    title(titulo, 'FontSize', 14)
    hold on
    y_fit=polyval(current_poly, x_fit);
    plot(x_fit, y_fit, '-g', 'LineWidth', 2)
    if idx_value==1
        eq_poly = sprintf('y = %4fx + %.4f', current_poly(1), current_poly(2));
    elseif idx_value==2
        eq_poly = sprintf('y = %.4fx^2 + %4fx + %.4f', current_poly(1), current_poly(2), current_poly(3));
    elseif idx_value==3
        eq_poly = sprintf('y = %.4fx^3 + %.4fx^2 + %4fx + %.4f', current_poly(1), current_poly(2), current_poly(3), current_poly(4));
    end
    text_x = min(x_fit) + 0.1*(max(x_fit) - min(x_fit));
    text_y = min(y_fit) + 0.05*(max(y_fit) - min(y_fit));
    text(text_x, text_y, eq_poly, 'FontSize', 12, 'BackgroundColor', 'w')
    if m==1
        R2_text=sprintf('R²= %.4f', R2_540);
    elseif m==2
        R2_text=sprintf('R²= %.4f', R2_595);
    end
    text_x2 = max(x_fit) - 0.15*(max(x_fit) - min(x_fit));
    text_y2 = min(y_fit) + 0.15*(max(y_fit) - min(y_fit));
    text(text_x2, text_y2, R2_text, 'FontSize', 12, 'BackgroundColor', 'w')
    hold off
    filename=['output/bradford_standard_curve_' wavelength];
    print ('-dpng', filename)
end

% calculating samples concentrations and matching it into the curve -
% Newton-Raphson
final.dc540_roots=[];
final.di540_roots=[];
final.dc595_roots=[];
final.di595_roots=[];
for n=1:2
    if n==1
        wavelength='540';
    else
        wavelength='595';
    end
    current_coef=['coef_' wavelength];
    current_poly=value.(current_coef);
    current_linear_coef=['linear_coef_' wavelength];
    current_lin=lin.(current_linear_coef);
    for k=1:2
        if k==1
            sample_type='c';
        else
            sample_type='i';
        end
        current_sample=['d' sample_type wavelength '_mean'];
        current_abs=sample.(current_sample);
        current_group=['d' sample_type wavelength '_roots'];
        for l=1:length(current_abs)
            error_absol=1000;
            f=current_poly;
            disp('f')
            disp(roots(f))
            f(end)=f(end)-current_abs(l);
            disp('root')
            disp(['d' sample_type wavelength '_' string(l)])
            disp(roots(f))
            df=polyder(f);
            current_lin_root=(current_lin(2)-current_abs(l))/current_lin(1);
            x=current_lin_root;
            counter=0;
            while error_absol>1e-4 && counter<100
                f_x=polyval(f, x);
                df_x=polyval(df,x);
                xi=x-(f_x/df_x);
                error_absol=abs(xi-x);
                x=xi;
                counter=counter+1;
                if x<0 && current_lin_root>0
                    while x<0
                        x=rand()*current_lin_root;
                    end
                end
            end
            if counter<100
                final.(current_group)(end+1)=xi;
            else
                final.(current_group)(end+1)=current_lin_root;
                disp('iteration out')
                disp(final.(current_group)(end))
            end
        end
    end
end

disp('dc 540nm roots:')
disp(final.dc540_roots)
disp('di 540nm roots:')
disp(final.di540_roots)
disp('relative difference between dc540 and di540 (%):')
relative_540=abs((final.dc540_roots-final.di540_roots)./final.dc540_roots)*100;
disp(relative_540)
disp('dc 595nm roots:')
disp(final.dc595_roots)
disp('di 595nm roots:')
disp(final.di595_roots)
disp('relative difference between dc595 and di595 (%):')
relative_595=abs((final.dc595_roots-final.di595_roots)./final.dc595_roots)*100;
disp(relative_595)

% saving manipulated data
result_roots=[final.dc540_roots; final.di540_roots; relative_540; final.dc595_roots; final.di595_roots; relative_595];
file_name=['output/sample_concentration.csv'];
fid = fopen(file_name, 'w');
fprintf(fid,'sample concentration data; \n ');
fprintf(fid,' \n');
fprintf(fid,'dc 540 (mg/mL);;; di 540 (mg/mL);;; relative difference dc/di (percentage);;; dc 595 (mg/mL);;; di 595 (mg/mL);;; relative difference dc/di (percentage);;; \n');
fprintf(fid,'%6.4f;', result_roots');
fclose(fid);

% plotting target concentrations on the original graph

for m=1:2
    if m==1
        wavelength='540';
    else
        wavelength='595';
    end
    for n=1:3
        if n==1
            sample_type='c';
        else
            sample_type='i';
        end
        current_coef=['coef_' wavelength];
        current_poly=value.(current_coef);
        current_white=['white' wavelength '_mean'];
        current_data=data.(current_white);
        current_sample=['d' sample_type wavelength '_mean'];
        current_abs=sample.(current_sample);
        current_group=['d' sample_type wavelength '_roots'];
        current_roots=final.(current_group);
        figure
        plot(albu_standard, current_data, 'o', 'LineWidth', 2)
        grid on
        xlabel('mg/mL', 'FontSize', 12)
        ylabel('abs', 'FontSize', 12)
        max_white_mean=max(current_data);
        ylim([0 max_white_mean+0.01])
        titulo=['Bradford curve with samples - ' wavelength 'nm'];
        title(titulo, 'FontSize', 14)
        hold on
        y_fit=polyval(current_poly, x_fit);
        plot(x_fit, y_fit, '-g', 'LineWidth', 2)
        handle1=plot(current_roots, current_abs, 's', 'MarkerEdgeColor','k','MarkerFaceColor','r', 'MarkerSize', 8);
        if n==3
            second_sample=['dc' wavelength '_mean'];
            second_abs=sample.(second_sample);
            second_group=['dc' wavelength '_roots'];
            second_roots=final.(second_group);
            handle2=plot(second_roots, second_abs, 's', 'MarkerEdgeColor','k','MarkerFaceColor','b', 'MarkerSize', 8);
            current_roots=[current_roots second_roots];
            current_abs=[current_abs second_abs];
        end
        for o=1:length(current_roots)
            x=current_roots(o);
            y=current_abs(o);
            line([x x], [0 y], 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
            line([0 x], [y y], 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
            if n==3 && o>3
                x2=current_roots(o);
                y2=current_abs(o);
                line([x2 x2], [0 y2], 'Color', [1 0 1], 'LineStyle', '--')
                line([0 x2], [y2 y2], 'Color', [1 0 1], 'LineStyle', '--')
            end
        end
        if n==3
            legend([handle1, handle2], {'di', 'dc'}, 'Location', 'southeast');
        else
            legend(handle1, ['d' sample_type], 'Location', 'southeast');
        end
        min_root=min(current_roots);
        max_root=max(current_roots);
        if min_root>0
            if max_root>max_albu
                xlim([0 max_root+0.005*max_root])
            else
                xlim([0 max_albu+0.005*max_albu])
            end
        else
            if max_root>max_albu
                xlim([min_root-0.005*abs(min_root) max_root+0.005*max_root])
            else
                if abs(max_albu-max_root)<0.15*max_root
                    xlim([min_root-0.005*abs(min_root) max_albu+0.005*max_albu])
                else
                    % change % position to generate a zoomed picture or
                    % standard version
                    xlim([min_root-0.005*abs(min_root) max_albu+0.005*max_albu])
                    %xlim([min_root-0.005*abs(min_root) max_root+0.1*max_root])
                end
            end
        end
        hold off
        filename=['output/bradford_curve_sample_' wavelength '-' num2str(n)];
        print ('-dpng', filename)
    end
end

% statistical analysis

% pearson correlation for dilution vs absorbance
sample_volumes=[100 50 20];
sample_volumes=sample_volumes*1e-6;
water_volume=[0 50 80];
water_volume=water_volume*1e-6;
bradford_volume=1e-3;
dilution_volume=sample_volumes ./ (sample_volumes+water_volume+bradford_volume);
[r1, p1]=corrcoef(sample.dc540_mean, dilution_volume);
[r2, p2]=corrcoef(sample.di540_mean, dilution_volume);
[r3, p3]=corrcoef(sample.dc595_mean, dilution_volume);
[r4, p4]=corrcoef(sample.di595_mean, dilution_volume);
disp('dc 540nm dilution correlation:')
disp(r1(1,2))
disp(p1(1,2))
disp('di 540nm dilution correlation:')
disp(r2(1,2))
disp(p2(1,2))
disp('dc 595nm dilution correlation:')
disp(r3(1,2))
disp(p3(1,2))
disp('di 595nm dilution correlation:')
disp(r4(1,2))
disp(p4(1,2))

% t-test for controlled/induced comparison and wavelength relevance

[h1, p5]=ttest(sample.dc540_mean, sample.di540_mean);
[h2, p6]=ttest(sample.dc595_mean, sample.di595_mean);
[h3, p7]=ttest2(sample.dc540_mean, sample.dc595_mean);
[h4, p8]=ttest2(sample.di540_mean, sample.di595_mean);
disp('dc/di 540nm correlation t:')
disp(h1)
disp(p5)
disp('dc/di 595 correlation t:')
disp(h2)
disp(p6)
disp('dc 540nm/595nm correlation t2:')
disp(h3)
disp(p7)
disp('di 540nm/595nm correlation t2:')
disp(h4)
disp(p8)

% white standard curve linearity and wavelength test

[h5, p9]=ttest(data.white540_mean, data.white595_mean);
[r5, p10]=corrcoef(albu_standard, data.white540_mean);
[r6, p11]=corrcoef(albu_standard, data.white595_mean);
disp('white standard wavelength correlation:')
disp(h5)
disp(p9)
disp('white 540nm linearity strength:')
disp(r5(1,2))
disp(p10(1,2))
disp('white 595nm linearity strength:')
disp(r6(1,2))
disp(p11(1,2))

% saving statistics data
result_stats=[r1(1,2); r2(1,2); r3(1,2); r4(1,2); h1; h2; h3; h4; h5; r5(1,2); r6(1,2)];
p_results=[p1(1,2); p2(1,2); p3(1,2); p4(1,2); p5; p6; p7; p8; p9; p10(1,2); p11(1,2)];
file_name=['output/statistical_curves_analysis.csv'];
fdd = fopen(file_name, 'w');
fprintf(fdd,'experimental curves statistical data; \n ');
fprintf(fdd,' \n');
fprintf(fdd,['data type; dc/540 pearson dilution; di/540 pearson dilution; dc/595 pearson dilution; ' ...
    'di/595 pearson dilution; dc/di 540 t; dc/di 595 t; dc 540/595 t2; di 540/595 t2; white t; white 540 pearson; ' ...
    'white 595 pearson; \n']);
fprintf(fdd,'results; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; \n', result_stats);
fprintf(fdd,'pvalues; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f; %6.4f;', p_results);
fclose(fdd);
