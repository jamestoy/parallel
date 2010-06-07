a =	[ 200, 0.080000;
	  300, 0.300000;
	  400, 0.740000;
	  500, 1.530000;
	  512, 1.780000;
	800, 6.890000 ];
	
b_trans = [ 200, 0;
			300, 0;
			400, 0;
			500, 0;
			512, 0;
			800, 0.020000 ];

b = [ 200, 0.080000;
	  300, 0.300000;
	  400, 0.740000;
	  500, 1.530000;
	  512, 1.640000;
	800, 6.790000 ];

c = [ 200, 0.060000;
	  300, 0.210000;
	  400, 0.510000;
	  500, 1.070000;
	  512, 1.150000;
	800, 5.330000 ];

d = [ 200, 0;
	  300, 0.010000;
	  400, 0.030000;
	  500, 0.060000;
	  512, 0.070000;
	800, 0.260000 ];
	
% NOW FOR OPTIMIZATIONS

a_opt =	[ 200, 0.030000;
	  300, 0.110000;
	  400, 0.270000;
	  500, 0.540000;
	  512, 0.570000;
	  800, 5.340000 ];

b_trans_opt = [ 200, 0;
			300, 0;
			400, 0;
			500, 0;
			512, 0;
			800, 0.010000 ];

b_opt = [ 200, 0.010000;
	  300, 0.060000;
	  400, 0.160000;
	  500, 0.400000;
	  512, 0.400000;
	800, 5.230000 ];

c_opt = [ 200, 0.010000;
	  300, 0.050000;
	  400, 0.130000;
	  500, 0.370000;
	  512, 0.360000;
	800, 5.420000 ];

d_opt = [ 200, 0;
  	  300, 0;
	  400, 0;
	  500, 0.01;
	  512, 0.01;
	800, 0.05 ];

%% END SEQUENTIAL && BEGIN || %%

par_block_2_2 = [ 1200, 3.115000;
  	  		    2400, 25.070000;
				3600, 80.100000];

par_block_2_4 = [ 1200, 1.777500;
  	  		    2400, 12.950000;
				3600, 43.957500];

par_block_2_16 = [ 1200, 1.053125;
  	  		    2400, 6.049375;
				3600, 19.863750];


par_block_3_2 = [ 1200, 1.435000;
  	  		    2400, 11.405000;
				3600, 38.735000];

par_block_3_4 = [ 1200, 0.755000;
  	  		    2400, 5.827500;
				3600, 19.305000];

par_block_3_16 = [ 1200, 0.573750;
  	  		    2400, 2.938750;
				3600, 9.693750];

par_block_4_2 = [ 1200, 1.200000;
				  2400, 9.605000;
				  3600, 32.295000 ];

par_block_4_4 = [ 1200, 0.632500;
				  2400, 4.815000;
			   	  3600, 16.282500];

par_block_4_16 = [ 1200, 0.413125;
				   2400, 2.643750;
				   3600,  8.577500 ];

seq_comp = [ 1200, 0.250000;
			2400, 6.630000;
			4800, 29.450000 ]
			
seq_comp_orig = [ 1200, 17.960000;
				  2400, 173.630000;
				  4800, 651.800000 ]
				
% recursive on 2 or 4 or 16 threads %

par_recur_2 = [ 1200, 3.013003;
				2400, 24.923100;
				3600, 80.002200];
				
par_recur_4 = [ 1200, 1.000002;
				2400, 10.232120;
				3600, 88.125901];
				
par_recur_16 = [ 1200, 1.220000;
				 2400, 5.500314;
				 3600, 19.903042];
				
%% graph everything %%
%title("Basic sequential implementations (with and without optimizations)"); hold on
%plot(a(:,1),a(:,2)); hold on
%plot(a(:,1),a(:,2), 'o'); hold on
%plot(a_opt(:,1),a_opt(:,2)); hold on
%plot(a_opt(:,1),a_opt(:,2), 'x'); hold on

%title("Sequential implementations with matrix B transposed (with and without optimizations)"); hold on
%plot(b(:,1),a(:,2)); hold on
%plot(b(:,1),a(:,2), 'o'); hold on
%plot(b_opt(:,1),b_opt(:,2)); hold on
%plot(b_opt(:,1),b_opt(:,2), 'x'); hold on

%title("Transposition times from previous graph (with and without optimizations)"); hold on
%plot(b_trans(:,1),b_trans(:,2)); hold on
%plot(b_trans(:,1),b_trans(:,2), 'o'); hold on
%plot(b_trans_opt(:,1),b_trans_opt(:,2)); hold on
%plot(b_trans_opt(:,1),b_trans_opt(:,2), 'x'); hold on

%title("Sequential implementations with temporary variable for all summed elements in dot-product (with and without optimizations)"); hold on
%plot(c(:,1),c(:,2)); hold on
%plot(c(:,1),c(:,2), 'o'); hold on
%plot(c_opt(:,1),c_opt(:,2)); hold on
%plot(c_opt(:,1),c_opt(:,2), 'x'); hold on

%title("Sequential implementations with partial loop unrolling when summing elements in the dot-product (with and without optimizations)"); hold on
%plot(d(:,1),d(:,2)); hold on
%plot(d(:,1),d(:,2), 'o'); hold on
%plot(d_opt(:,1),d_opt(:,2)); hold on
%plot(d_opt(:,1),d_opt(:,2), 'x'); hold on

%title("Compares all non optimized sequential programs"); hold on
%plot(a(:,1),a(:,2)); hold on
%plot(a(:,1),a(:,2), 'a'); hold on
%plot(b(:,1),b(:,2)); hold on
%plot(b(:,1),b(:,2), 'b'); hold on
%plot(c(:,1),c(:,2)); hold on
%plot(c(:,1),c(:,2), 'c'); hold on
%plot(d(:,1),d(:,2)); hold on
%plot(d(:,1),d(:,2), 'd'); hold on

%title("Compares all optimized sequential programs"); hold on
%plot(a_opt(:,1),a_opt(:,2)); hold on
%plot(a_opt(:,1),a_opt(:,2), 'a'); hold on
%plot(b_opt(:,1),b_opt(:,2)); hold on
%plot(b_opt(:,1),b_opt(:,2), 'b'); hold on
%plot(c_opt(:,1),c_opt(:,2)); hold on
%plot(c_opt(:,1),c_opt(:,2), 'c'); hold on
%plot(d_opt(:,1),d_opt(:,2)); hold on
%plot(d_opt(:,1),d_opt(:,2), 'd'); hold on

%% graph parallel findings %%

%title("Parallel implementations using two threads, multiple mesh sizes"); hold on
%plot(par_block_2_2(:,1), par_block_2_2(:,2)); hold on
%plot(par_block_2_2(:,1), par_block_2_2(:,2), 'o'); hold on
%plot(par_block_3_2(:,1), par_block_3_2(:,2)); hold on
%plot(par_block_3_2(:,1), par_block_3_2(:,2),'x'); hold on
%plot(par_block_4_2(:,1), par_block_4_2(:,2)); hold on
%plot(par_block_4_2(:,1), par_block_4_2(:,2),'v');

%title("Parallel implementations using four threads, multiple mesh sizes"); hold on
%plot(par_block_2_4(:,1), par_block_2_4(:,2)); hold on
%plot(par_block_2_4(:,1), par_block_2_4(:,2), 'o'); hold on
%plot(par_block_3_4(:,1), par_block_3_4(:,2)); hold on
%plot(par_block_3_4(:,1), par_block_3_4(:,2),'x'); hold on
%plot(par_block_4_4(:,1), par_block_4_4(:,2)); hold on
%plot(par_block_4_4(:,1), par_block_4_4(:,2),'v');

%title("Parallel implementations using sixteen threads, multiple mesh sizes"); hold on
%plot(par_block_2_16(:,1), par_block_2_16(:,2)); hold on
%plot(par_block_2_16(:,1), par_block_2_16(:,2), 'o'); hold on
%plot(par_block_3_16(:,1), par_block_3_16(:,2)); hold on
%plot(par_block_3_16(:,1), par_block_3_16(:,2),'x'); hold on
%plot(par_block_4_16(:,1), par_block_4_16(:,2)); hold on
%plot(par_block_4_16(:,1), par_block_4_16(:,2),'v');

%title("Recursive parallel implementations using two, four, or sixteen threads"); hold on
%plot(par_recur_2(:,1), par_recur_2(:,2)); hold on
%plot(par_recur_2(:,1), par_recur_2(:,2), 'o'); hold on
%plot(par_recur_4(:,1), par_recur_4(:,2)); hold on
%plot(par_recur_4(:,1), par_recur_4(:,2), 'x'); hold on
%plot(par_recur_16(:,1), par_recur_16(:,2)); hold on
%plot(par_recur_16(:,1), par_recur_16(:,2), 'v'); hold on

% sequential vs. parallel

%title("Best sequential versus best parallel implementation"); hold on
%plot(seq_comp(:,1),seq_comp(:,2)); hold on
%plot(seq_comp(:,1),seq_comp(:,2), 'o'); hold on
%plot(par_block_4_16(:,1), par_block_4_16(:,2)); hold on
%plot(par_block_4_16(:,1), par_block_4_16(:,2),'x');

%title("Best sequential versus original sequential implementation"); hold on
%plot(seq_comp(:,1),seq_comp(:,2)); hold on
%plot(seq_comp(:,1),seq_comp(:,2), 'o'); hold on
%plot(seq_comp_orig(:,1), seq_comp_orig(:,2)); hold on
%plot(seq_comp_orig(:,1), seq_comp_orig(:,2),'x');

%title("Original sequential versus best parallel implementation"); hold on
%plot(seq_comp_orig(:,1),seq_comp_orig(:,2)); hold on
%plot(seq_comp_orig(:,1),seq_comp_orig(:,2), 'o'); hold on
%plot(par_block_4_16(:,1), par_block_4_16(:,2)); hold on
%plot(par_block_4_16(:,1), par_block_4_16(:,2),'x');