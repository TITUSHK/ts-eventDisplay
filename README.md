# ts-eventDisplay
Repository for event display

The program titus_eve_4Reco.C will draw events pruduced by ts-WChRecoSandBox.
To run it,start root and then do
.x titus_eve_4Reco.C++

OR

.L titus_eve_4Reco.C++
and then
titus_eve_4Reco()

Next time you run root do

.L titus_eve_4Reco.so
and then do
titus_eve_4Reco()

At the first dialog choose your file of ts-WChRecoSandBox events, at the second 
you can choose the matching file of reconstruction output.

Alternatively, specify the input file on the command line, e.g.:

titus_eve_4Reco("nu_antinue_1066_out_12in.root","nu_antinue_1066_12in_reco_4.root")

There is a 3D display and an unrolled display.
You can select the event.
Look the Eve tab to do things like selecting what is drawn, and changing
the cuts on the hits. 

Hits are drawn with a colour to represent time and the size of the circle to
represent the total charge seen on a given PMT.The hits in the 3D event are under Event, those in the
unrolled view are under Scenes>Unrolled Event>PMT Hits(Unrolled). Once hits are 
selected you can control them using the Style tab, and histogram the value (time)
using the Info tab. 
In 3D neutron captures are shown separately as a yellow cube.

To reduce the number of truth particles shown choose
Truth Particles, click on top level objects and turn off the tick next to "Children".

There is also a global control to just show primary tracks, or not.

Reconstructed tracks from the "Final" sections are shown as straight lines.
Objects with no direction just as cubes at the vertex.

Contact Alex Finch (A.Finch@lancaste.ac.uk) for bug reports and update requests.






