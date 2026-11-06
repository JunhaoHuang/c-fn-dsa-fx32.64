#define main full_benchmark_main
#include "../speed_masked.c"
#undef main

int
main(void)
{
	use_tsc = 1;
	counter_init();
	gadget_ctx context;
	init_gadget_ctx(&context);
	print_mod_group(&context, 12289);
	print_mod_group(&context, 18433);
	print_mod_group(&context, UINT32_C(226523137));
	print_mq_backend_group(9);
	print_mq_backend_group(10);
	counter_close();
	return 0;
}
