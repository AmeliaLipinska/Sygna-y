import sys
import os

sys.path.insert(0, os.path.abspath("build/Release"))
import signals

def main():
    valid_types = {"sin", "cos", "square", "sawtooth"}
    while True:
        type1 = input("Typ sygnalu 1 (sin/cos/square/sawtooth): ").lower()
        if type1 in valid_types:
            break
        print("Bledny typ, wybierz ponownie.")

    while True:
        try:
            freq1 = float(input("Czestotliwosc 1 (Hz): "))            
            amp1 = float(input("Amplituda 1: "))
            break
        except ValueError:
            print("Zla wartosc, wpisz ponownie.")

    while True:
        type2 = input("Typ sygnalu 2 (sin/cos/square/sawtooth): ").lower()
        if type2 in valid_types:
            break
        print("Bledny typ, wybierz ponownie.")

    while True:
        try:
            freq2 = float(input("Czestotliwosc 2 (Hz): "))            
            amp2 = float(input("Amplituda 2: "))
            break
        except ValueError:
            print("Zla wartosc, wpisz ponownie.")

    while True:
        try:
            t_start = float(input("Podaj czas startowy: "))
            t_end = float(input("Podaj czas koncowy: "))
            num_samples = int(input("Podaj ilosc probek: "))
            threshold = float(input("Podaj prog: "))

            break
        except ValueError:
            print("Zla wartosc, wpisz ponownie.")


    signals.run(type1, freq1, amp1, type2, freq2, amp2, t_start, t_end, num_samples, threshold)


if __name__ == "__main__":
    main()