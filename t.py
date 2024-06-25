
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', default="0" , type=str)
    parser.add_argument('--coarse', default=0 , type=float)
    args = parser.parse_args()

    

    print(args.id)


