import os, glob, sys, math

methods = ["dirichlet", "asap", "symdirichlet", "arap"]
start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
step = int(sys.argv[2]) if len(sys.argv) > 2 else 1
end = int(sys.argv[3]) if len(sys.argv) > 3 else None
num = math.floor((len(glob.glob("test_data/*.off")[:end])-1-start)/step)

for k,i in enumerate(sorted(glob.glob("test_data/*.off"), key = os.path.getsize)[start:end:step]):
    print(k,"of", num)
    for m in methods:
        print(f"build/cli {m} {i} output")
        try:
            sts = os.system(f"build/cli {m} {i} output")
            if(sts == 2):
                sys.exit(0)
        except Exception as e:
            print("Failed method",m,"for mesh",i)
        print()
    print()
    print()
