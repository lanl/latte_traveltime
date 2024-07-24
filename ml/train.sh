
# Train 2D inferrence NN
python main2_infer.py \
	--ntrain=2000 --nvalid=100 \
	--n1=256 --n2=256 \
	--batch_train=8 --batch_valid=8 --lr=0.5e-4 \
	--gpus_per_node=1 --epochs=50 \
	--dir_output=result2_infer 

# Train 2D refinement NN
python main2_refine.py \
	--ntrain=2000 --nvalid=100 \
	--n1=256 --n2=256 \
	--batch_train=8 --batch_valid=8 --lr=0.5e-4 \
	--gpus_per_node=1 --epochs=50 \
	--dir_output=result2_refine 

# Train 3D inferrence NN
python main3_infer.py \
	--ntrain=1000 --nvalid=100 \
	--n1=128 --n2=180 --n3=180 \
	--batch_train=1 --batch_valid=1 --lr=0.5e-4 \
	--gpus_per_node=4 --nodes=2 \
	--dir_output=result3_infer

# Train 3D refinement NN
python main3_refine.py \
	--ntrain=1000 --nvalid=100 \
	--n1=128 --n2=180 --n3=180 \
	--batch_train=1 --batch_valid=1 --lr=0.5e-4 \
	--gpus_per_node=4 --nodes=2 \
	--dir_output=result3_refine
