//Code provided by Rachelle-Heim Boissier and Yann Rotella
//Everything is free of use, re-use and change
//Don't hesitate to contact us if you find flaws
//Or more importantly if you find ways to improve the attack !!!

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int b = 100; //keccak variant of width 100
int nr = 2; //2 round attack
int r = 30; //rate (generic attack in 2^30)
int c = 70; //capacity 
int w = 4;  //the omega of keccak
int l = 2;  //the l of keccak

//initializes the state to the all 0 value
void init_zero(int state[5][5]){
	for (int x = 0; x < 5; ++x)
	{
		for (int y = 0; y < 5; ++y)
		{
			state[x][y] = 0;
		}
	}
}

//rotate an int of length w
int rotate(int n, int d){
	return ((n >> d)|(n << (w-d)) & 0xf);
}

//theta mapping
void theta(int state[5][5]){
	int c[5] = {0x0,0x0,0x0,0x0,0x0};
	for (int x = 0; x < 5; ++x)
	{
		for (int y = 0; y < 5; ++y)
		{
			c[x] ^= state[x][y];
		}
	}
	for (int x = 0; x < 5; ++x)
	{
		for (int y = 0; y < 5; ++y)
		{
			state[x][y] ^= (c[(x+4)%5]^(rotate(c[(x+1)%5],3)));
		}
	}
}

//rho mapping
void rho(int state[5][5]){
	int x = 1;
	int y = 0;
	for (int t = 0; t < 24; ++t)
	{
		int rot = ((t+1)*(t+2)/2)%w;
		state[x][y] = rotate(state[x][y],w-rot);
		int tmp = y;
		y = (2*x + 3*tmp)%5;
		x = tmp;
	}
}

//pi mapping
void pi(int state[5][5]){
	int tmp[5][5];
	for (int x = 0; x < 5; ++x)
	{
		for (int y = 0; y < 5; ++y)
		{
			tmp[x][y] = state[(x+3*y)%5][x];
		}
	}
	for (int x = 0; x < 5; ++x)
	{
		for (int y = 0; y < 5; ++y)
		{
			state[x][y] = tmp[x][y];
		}
	}
}

//chi mapping
void chi(int state[5][5]){
	for (int y = 0; y < 5; ++y)
	{
		int tmp[5];
		for (int x = 0; x < 5; ++x)
		{
			tmp[x] = state[x][y] ^ ((0xf ^ state[(x+1)%5][y]) & (state[(x+2)%5][y]));
		}
		for (int x = 0; x < 5; ++x)
		{
			state[x][y] = tmp[x];
		}
	}
}

//add constant
void iota(int state[5][5], int ir){
	if (ir == 14)
	{
		state[0][0] ^= 0x9;
	}
	else
	{
		state[0][0] ^= 0x3;
	}		
}

//messages are stored in 30 bits bocks (rate) while state parts are stored in 4 bits
void conversion(int msg, int state[5][5]){
	int tmp = msg;
	int i = 0;
	while (tmp > 0)
	{
		state[(i/w)%5][i/(w*5)] ^= (tmp & 0xf);
		tmp = tmp >> 4;
		i += 4;
	}
}

//msg is a tab of 30 bits message blocks, absorb every 2 rounds
void absorb(int* msg, int length, int nr, int state[5][5]){
	for (int i = 0; i < length; ++i)
	{
		conversion(msg[i],state);
		for (int ir = (12 + 2*l - nr); ir < (12 + 2*l); ++ir)
		{
			theta(state);
			rho(state);
			pi(state);
			chi(state);
			iota(state, ir);
		}
	}
}

//Linear system derived for the attack (see Python programm)
int LinSys[33][3] = {{0x2a,0x04175471,0xe0e0a625}, {0x30,0x26080242,0x20242b26}, {0x3,0x04138041,0x890090b0}, {0x3,0x46235530,0x81d68431},
{0x23,0x32530257,0x2180a717}, {0x3a,0x3482823b,0xa8a647a2}, {0x32,0x6700de22,0x28e18382}, {0x0,0x22,0x00002000},
{0x2,0x04031150,0x80b48421}, {0x39,0x780106e3,0x21ee4712}, {0x3,0x56035410,0x81c68430}, {0x2a,0x34038210,0xa884e6a2},
{0x3,0x06031130,0x81968435}, {0x21,0xd2004d00,0x21f22219}, {0x2,0x04030110,0x80b48421}, {0x1,0x0401a040,0x8b008090},
{0x7,0x0403c050,0x8d8080f0}, {0xb,0x4603c510,0x89d684b1}, {0x2a,0x34038232,0xb8a6c6a2}, {0x2,0x04011110,0x80b48421},
{0x0,0x1010d461,0x08440485}, {0x0,0x4202c510,0x08d60481}, {0x0,0x0,0x0}, {0x0,0x0,0x0},
{0x20,0x36128253,0xa88067a7}, {0x0,0x0,0x0}, {0x10,0x6203cf22,0x885425a3}, {0x3,0x06031112,0x8196a431},
{0x23,0x02031130,0xa1928635}, {0x3a,0x74028633,0xa8e6c7a2}, {0x2a,0x0417c531,0xe8f4a2a4}, {0x0,0x1000,0x00000010},
{0x0,0x20000,0x20000000}};

//matrix vector multiplication to give messages for building collision
int eval(int state[5][5],int poly[3]){
	int img = 0;
	for (int j = 0; j < 70; ++j)
	{
		int pos = j + r;
		img ^= ((poly[2 - (j/32)] >> (j%32)) & 1) & ((state[(pos/w)%5][pos/(w*5)]) >> (pos%w));
	}
	return img;
}

//from a state, give next message block to absorb and lie into subset X
int derive_message(int state[5][5]){
	int msg = 0;
	for (int i = 0; i < r; ++i)
	{
		msg ^= (eval(state,LinSys[i]) << i);
	}
	return msg;
}

//equivalent messages possible as the matrix is not full rank
void add_kernel(int e,int state[5][5]){
	state[0][1] ^= (e & 1) << 2;
    state[0][1] ^= ((e >> 1) & 1) << 3;
    state[1][1] ^= ((e >> 2) & 1) << 1;
    state[0][0] ^= (e & 1) << 2;
    state[0][0] ^= ((e >> 1) & 1) << 3;
    state[1][0] ^= ((e >> 2) & 1) << 1;
}

//test the probability after the chi mapping
int probabilistic_after_chi(int state[5][5]){
	int b0 = 1 ^(state[0][1] & 1)^(state[0][3] & 1);
    int b1 = 1 ^(state[0][3] & 1)^(state[0][4] & 1);
    int b2 = ((1<<2) ^ (state[2][0] & (1<<2)) ^ (state[2][1] & (1<<2))) >> 2;
    int b3 = ((1<<2) ^ (state[2][1] & (1<<2)) ^ (state[2][3] & (1<<2))) >> 2;
    int b4 = ((1<<3) ^ (state[3][1] & (1<<3)) ^ (state[3][2] & (1<<3))) >> 3;
    int b5 = ((1<<3) ^ (state[3][2] & (1<<3)) ^ (state[3][4] & (1<<3))) >> 3;
    return (b0 & b1 & b2 & b3 & b4 & b5);
}

//store capacity in 70 bit integers until you find a collision
//now it is on 60 bits to have something that can be run "almost fast"
__uint128_t capacity(int state[5][5]){
	int rot = 0;
	__uint128_t tmp = 0;
	int pos = 40;
	for (int i = 0; i < 15; ++i)
	{
		tmp ^= (__uint128_t) state[(pos/w)%5][pos/(w*5)] << (rot);
		rot += 4;
		pos += 4;
	}
	//tmp ^= (__uint128_t) (state[(30/w)%5][30/(w*5)] & 0xc) << (rot);
	//if you want full 70 bit capacity
	return tmp;
}

__uint64_t capacity64(int state[5][5]){
	int rot = 0;
	__uint64_t tmp = 0;
	int pos = 40;
	for (int i = 0; i < 10; ++i)
	{
		tmp ^= (__uint64_t) (state[(pos/w)%5][pos/(w*5)] << (rot));
		rot += 4;
		pos += 4;
	}
	return tmp;
}

//deep copy the state
void my_deep_copy(int state[5][5], int tab[5][5]){
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			tab[i][j] = state[i][j];
		}
	}
}

//used for quick-sort if you mount attack with memory
int compare(const void* a, const void* b)
{
    if(*(__uint64_t*)a - *(__uint64_t*)b < 0)
        return -1;
    if(*(__uint64_t*)a - *(__uint64_t*)b == 0)
        return 0;
    if(*(__uint64_t*)a - *(__uint64_t*)b > 0)
        return 1;
}

//the following is for generating pseudo-random messages
static inline int rotl(const int x, int k) {
	return (x << k) | (x >> (32 - k));
}
static int s[4];
int next(void) {
	const int result = s[0] + s[3];
	const int t = s[1] << 9;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = rotl(s[3], 11);
	return result;
}

int MESSAGE[5] = {0,0,0,0,0};

void cap_to_message(__uint128_t c){
	//initiate the prng with capacity information
	s[0] = (int) (c);
	s[1] = (int) (c>>32);
	s[2] = (int) (c>>64);
	s[3] = 0;
	for (int i = 0; i < 5; ++i)
		{
			MESSAGE[i] = (next() & 0x3fffffff);
		}
}

__uint128_t next_capacity(void) {
		int state[5][5];
		init_zero(state);
		absorb(MESSAGE,5,nr,state);
		int test = eval(state,LinSys[30]) | eval(state,LinSys[31]) | eval(state,LinSys[32]);
		while (test&1){
			absorb(MESSAGE,1,nr,state);//this choice is arbitrary
			test = eval(state,LinSys[30]) | eval(state,LinSys[31]) | eval(state,LinSys[32]);
		}
		int new_block = derive_message(state);
		int rate = 0;
		for (int i = 0; i < r; ++i){
			rate ^= ((state[(i/w)%5][i/(w*5)] >> (i%w)) << i);
			if ((i%w == 0) && (i < 28)){
				state[(i/w)%5][i/(w*5)] = ((new_block >> i)&0xf);
			}
			if (i == 28){
				int v = state[(i/w)%5][i/(w*5)];
				state[(i/w)%5][i/(w*5)] = (((new_block >> i)&0xf)^(v&0xc));
			}
		}
		int true_message = new_block^rate;
		//apply two rounds
		theta(state);
		rho(state);
		pi(state);
		chi(state);
		iota(state,14);
		theta(state);
		rho(state);
		pi(state);
		chi(state);
		iota(state,15);
		__uint128_t r = capacity(state);
		return r;
}

__uint128_t my_f(__uint128_t c){
	cap_to_message(c);
	return next_capacity();
}

//memory-less attack
int floyd(void){
	srand(time(NULL));
	s[0] = rand();
	s[1] = rand();
	s[2] = rand();
	s[3] = rand();
	for (int i = 0; i < 5; ++i)
		{
			MESSAGE[i] = (next() & 0x3fffffff) ;
		}
	__uint128_t c = next_capacity();
	__uint128_t tortoise = my_f(c);
	__uint128_t hare = my_f(my_f(c));
	int n = 1;
	while (tortoise != hare){
		tortoise = my_f(tortoise);
		hare = my_f(my_f(hare));
		n += 3;
		if ((n%10000) == 0)
		{
			printf("%d\n", n/10000);
		}
	}
	printf("number of calls to f: %d\n",n);
	int mu = 0;
    tortoise = c;
    while (tortoise != hare){
        tortoise = my_f(tortoise);
        hare = my_f(hare);
        mu += 1;
        if (mu%10000 == 0)
		{
			printf("%d\n", mu/10000);
		}
    }
    int lam = 1;
    hare = my_f(tortoise);
    while (tortoise != hare){
        hare = my_f(hare);
        lam += 1;
        if (lam%10000 == 0)
		{
			printf("%d\n", lam/10000);
		}
    }
    printf("Starting point of the cycle: %d\n", mu);
    printf("Length of the period: %d\n",lam);
}

void main(int argc, char *argv[])
{
	printf("Experiments:\n");
	for (int i = 0; i < 5; ++i)
	{
		floyd();
	}
}
