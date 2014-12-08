using Estimation
using Gadfly
using DSP

#
# Create a test signal similar to Li & Stoica 96
#

# Li & Stoica fig4
spectral_lines_f = [0.04, 0.05, 0.08, 0.09, 0.195, 0.215]
spectral_lines_p = vec(repmat([pi/4], length(spectral_lines_f)))
spectral_lines_a = [10, 10, 10, 2, 4, 6]

# Simulation properties
N = 128
noise_mean = 0
noise_vari = 1
num_trials = 100
padd_times = 16

function create_example(N, spectral_lines_f, spectral_lines_p, spectral_lines_a)
    clean = vec(zeros(N, 1))
    n = 1.0:N
    for i in 1:length(spectral_lines_f)
        clean = clean + spectral_lines_a[i] * sin(2 * pi * spectral_lines_f[i] * n + spectral_lines_p[i])
    end
    return clean
end

clean_signal = create_example(N, spectral_lines_f, spectral_lines_p, spectral_lines_a)
noise = randn(N)

z = clean_signal + noise


#
# Create a plot of the true spectrum
#

# Plot true spectrum
yline = zeros(length(0.0:0.0001:0.5))
for i in 1:length(spectral_lines_f); yline[round(spectral_lines_f[i]*10000)] = spectral_lines_a[i]; end

# Save plot as a layer so you can overlay them all later
true_plot = layer(
    x=0.0:0.0001:0.5,
    y=yline,
    Geom.line(),
    Theme(default_color=color("black")))

# Plot the true spectrum
p0 = plot(true_plot,
    Guide.title("True spectrum   N = $N"),
    Guide.xlabel("Frequency (Hz)"),
    Guide.ylabel("Amplitude"),
    Scale.y_continuous(minvalue=0, maxvalue=12))


#
# Plot FFT of signal
#

p = periodogram(z, nfft=length(z)*padd_times)
fft_plot = layer(
        x=freq(p),
        y=sqrt(power(p)/(length(z)/2)),
        Theme(default_color=color("blue")),
        Geom.line);
p1 = plot(fft_plot, true_plot,
    Guide.title("FFT spectrum"),
    Guide.xlabel("Frequency (Hz)"),
    Guide.ylabel("Amplitude"),
    Scale.y_continuous(minvalue=0, maxvalue=12));


#
# Plot Kaiser window FFT
#

window = kaiser(N, 4/pi)
window_power = mean(window.^2)

p = periodogram(z, nfft=length(z)*padd_times, window=window)

kaiser_plot = layer(
    x=freq(p),
    y=sqrt(power(p)/(length(z)/2)/sqrt(window_power)),
    Theme(default_color=color("purple")),
    Geom.line)

p2 = plot(kaiser_plot, true_plot,
    Guide.title("FFT spectrum with Kaiser window"),
    Guide.xlabel("Frequency (Hz)"),
    Guide.ylabel("Amplitude"),
    Scale.y_continuous(minvalue=0, maxvalue=12))


#
# Plot APES method
#

f = 0.00:0.001:0.5
a = apes(z, f)

apes_plot = layer(
    x=f,
    y=abs(a),
Theme(default_color=color("red")),
    Geom.line)

p5 = plot(apes_plot,
    true_plot,
Guide.title("APES spectrum with M = $(int(N/2))"),
    Guide.xlabel("Frequency (Hz)"),
    Guide.ylabel("Amplitude"),
    Scale.y_continuous(minvalue=0, maxvalue=12))


#
# Connect plots and save
#

set_default_plot_size(24cm, 12cm)
Pt = Gadfly.gridstack([p0 p1; p2 p5])
# draw(PNG("Li-Stoica-Fig1.png", 28cm, 28cm), Pt)
#=draw(PDF("Li-Stoica-Fig4-Presentation.pdf", 26cm, 16cm), Pt)=#
draw(PNG("APES-example.png", 26cm, 16cm), Pt)
#=Pt=#
