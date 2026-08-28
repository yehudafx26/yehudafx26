function [fn, vb, xf, ff] = generateGestures(ges_num, Nf, SR)

re = 0.03/0.38;
mi = 0.055/0.38;
fa = 0.073/0.38;
sol = 0.093/0.38;

la = re;
si = mi;
doo = fa;
ree = sol;

res = re+mi/2;

%% open string and press string

if ges_num == 1
ff = 0.5*ones(Nf,1);   ff(1:Nf/2) = 0;
fn = 0.5 * sin(2*pi*4*(0:Nf-1)/SR);
vb = 0.2 * sin(2*pi*4*(0:Nf-1)/SR);
xf = 0.085*ones(Nf,1);
end

%% a scale

if ges_num == 2
sec = floor(Nf/8);
ff = 0.5*ones(Nf,1);
xf = zeros(Nf,1);
xf(1:2*sec) = 0.085;
xf(2*sec: 3*sec) = 0.13;
xf(3*sec: 4*sec) = 0.25;
xf(3*sec: 5*sec) = 0.25;
xf(5*sec: 6*sec) = 0.13;
xf(6*sec: Nf) = 0.085;
fn = 0.5*ones(Nf,1);
vb = 0.3 * sin(2*pi*4*(0:Nf-1)/SR);
end

%% finger up and down

if ges_num == 3

    sec = floor(Nf/8);
    unit = [ 0.3 * ones(sec, 1); zeros(sec, 1)];
    ff = repmat(unit, 4, 1);
    xf = re*ones(size(ff, 1),1);
    fn = 0.4 * ones(size(ff, 1),1);

    vunit = [-0.2 * ones(sec, 1); 0.2 * ones(sec, 1)];
    vb = repmat(vunit, 4, 1);

end

%% a Teochew scale

if ges_num == 4
    sec = floor(Nf/10);
    unit = [zeros(sec, 1); 0.3 * ones(sec, 1)];
    ff = repmat(unit, 5, 1);
    xf1 = linspace(re,re,sec); 
    xf2 = linspace(re,res,sec); 
    xf3 = linspace(res,re,sec); 
    xf4 = linspace(re,fa,sec); 
    xf5 = linspace(fa,fa,sec); 
    xf6 = linspace(fa,res,sec); 
    xf7 = linspace(res,re,sec); 
    xf8 = linspace(res,re,sec); 
    xf9 = linspace(re,fa,sec); 
    xf10 = linspace(fa,re,sec); 
    xf = [xf1,xf2,xf3,xf4,xf5,xf6,xf7,xf8,xf9,xf10];

    fn = 0.5 * ones(size(xf, 2),1);
    vb = 0.2 * ones(size(xf, 2),1);
end

%% slide + vibrato

if ges_num == 5

    ff = 0.4 * ones(Nf, 1);  
    fn = 0.5 * ones(Nf, 1);
    vbv = 0.2; 
    
    vib_freq = 6; vib_amp = 0.01;    
    xf_base = zeros(Nf, 1);
    v_mask = ones(Nf, 1); 
    vb = zeros(Nf, 1);    
    
    t = (0:Nf-1)' / SR; 
    p1 = floor(Nf/4);   
    p2 = floor(Nf/2);   
    slide_len = floor(p1 * 0.5); 
    
    acc_len = floor(p1 * 0.15); % acc time
    dec_len = floor(p1 * 0.15); % decay time

    gen_vb_smooth = @(total_len, v_max, dir, a_len, d_len) [ ...
        dir * v_max * sin(linspace(0, pi/2, a_len))'; ...           % 起弓加速 (0 -> max)
        dir * v_max * ones(total_len - a_len - d_len, 1); ...       % 稳定段 (max)
        dir * v_max * cos(linspace(0, pi/2, d_len))' ...            % 收弓减速 (max -> 0)
    ];

    % --- note 1 ---
    xf_base(1:p1) = re;
    vb(1:p1) = gen_vb_smooth(p1, vbv, 1, acc_len, dec_len); 

    % --- note 2 ---
    slide_idx1 = p1+1 : p1+slide_len;
    xf_base(slide_idx1) = linspace(re, fa, slide_len);
    v_mask(slide_idx1) = 0; 
    xf_base(p1+slide_len+1 : p2) = fa;
    vb(p1+1 : p2) = gen_vb_smooth(p2 - p1, vbv, -1, acc_len, dec_len);

    % --- note 3 ---
    slide_idx2 = p2+1 : p2+slide_len;
    xf_base(slide_idx2) = linspace(fa, sol, slide_len);
    v_mask(slide_idx2) = 0; 
    xf_base(p2+slide_len+1 : Nf) = sol;
    
    vb(p2+1 : Nf) = gen_vb_smooth(Nf - p2, vbv, 1, acc_len, dec_len);

    vibrato = (vib_amp * sin(2 * pi * vib_freq * t)) .* v_mask;
    xf = xf_base + vibrato;
    xf = max(xf, 0.01);

end

vb = vb(:);
fn = fn(:);
xf = xf(:);
ff = ff(:);

end
