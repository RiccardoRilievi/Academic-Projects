%% TABELLA TRASFERIMENTI

DeltaT_4 = deltat_4_CP + deltat_4_CPA +...
            deltat_4_P + Deltat_4_BT + deltat_4_f;
DeltaV_4 = abs(DeltaV_4_CP) + abs(DeltaV_4_CPA) + abs(DeltaV_4_BT);

figure(4)
plot(DeltaT_4, DeltaV_4,"o", LineWidth=5.5, color="#7E2F8E");
hold on;
grid on;
grid minor;

DeltaV_A3 = abs(DeltaV_tan) + abs(DeltaV_1_CP) + abs(DeltaV_CPA);
DeltaT_A3 = deltat_tan + deltat_1_CPA + deltat_2_CPA + deltat_1_CP;

figure(4)
plot(DeltaT_A3, DeltaV_A3, "o", LineWidth=5.5, color="#EDB120")

DeltaT_ALT1 = deltat_1_CPA + deltat_1_CP + deltat_sec;
DeltaV_ALT1 = abs(DeltaV_1_CPA) + abs(DeltaV_1_CP) + abs(DeltaV_sec);

figure(4)
plot(DeltaT_ALT1, DeltaV_ALT1,"o", LineWidth=5.5, color="#77AC30");

xlabel("Time (s)")
ylabel("Velocity (km/s)")
axis([0 35000 2 5])
legend("standard strategy", "first alternative strategy", "second alternative strategy")