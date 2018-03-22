ids_all = [];
ads_all = [];

for k = 1:length(ids_i)
   
    ids_all = [ids_all;F0PARAMS{ids_i(k)}];
    
end


for k = 1:length(ads_i)
   
    ads_all = [ads_all;F0PARAMS{ads_i(k)}];
    
end

figure(2);clf;
subplot(2,2,2);hold on;
scatter(ads_all(:,1),ads_all(:,2));
scatter(ids_all(:,1),ids_all(:,2))

legend({'ADS','IDS'})
xlim([-0.3 0.3]);
ylim([-2 2]);grid;
xlabel('polycoeff a_2');
ylabel('polycoeff a_1');
subplot(2,2,4);hold on;
mu = mean(ads_all(:,1));
sigma = std(ads_all(:,1));
x = -0.5:0.01:0.5;
y = normpdf(x,mu,sigma);

plot(x,y,'LineWidth',2);
mu = mean(ids_all(:,1));
sigma = std(ids_all(:,1));
x = -0.5:0.01:0.5;
y = normpdf(x,mu,sigma);
plot(x,y,'LineWidth',2);

xlim([-0.3 0.3]);
grid;
subplot(2,2,1);hold on;
mu = mean(ads_all(:,2));
sigma = std(ads_all(:,2));
x = -4:0.01:4;
y = normpdf(x,mu,sigma);

plot(x,y,'LineWidth',2);
mu = mean(ids_all(:,2));
sigma = std(ids_all(:,2));
x = -4:0.01:4;
y = normpdf(x,mu,sigma);
plot(x,y,'LineWidth',2);

xlim([-2 2]);
grid;

