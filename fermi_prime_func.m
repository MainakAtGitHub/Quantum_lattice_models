function f=fermi_prime_func(x)
        %Derivative (prime) of Fermi's function
f=-fermi_func(x).*fermi_func(-x);
    