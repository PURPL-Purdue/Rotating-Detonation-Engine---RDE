%example of the ceaRocket function being used
r = ceaRocket('ox','Air','fuel','H2','phi',1,'Pc',1000,'PcUnits','psia');
disp(r.isp); % m/s
disp(r.ivac); % m/s
disp(r.cstar); % m/s
disp(r.cf);
disp(r.AeAt);