Gas sensor array under dynamic gas mixtures - RAW coordinates
Source: UCI Machine Learning Repository, dataset 322, file ethylene_CO.txt
        (Fonollosa, 2015; DOI 10.24432/C5WP4C)

Columns: time_seconds  then  c0..c15 (16 RAW sensor channels, as recorded)
Each file = 2000 consecutive samples of a 100 Hz acquisition.

Time is the recorded clock, in seconds, and is NOT rescaled, shifted or otherwise
altered. As in the source it is quantised to 0.01 s, so some consecutive rows share
a value and others step by 0.02 s; that quantisation is the recording's own and is
left in place.

Channels are NOT rescaled: the ratio between them is as recorded.

Two columns present in the UCI file are dropped: the methane (or CO) and ethylene
concentration set points. They are experimental settings rather than sensor state,
so they are not part of the trajectory.
